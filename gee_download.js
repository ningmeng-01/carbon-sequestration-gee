var zhumadian = ee.FeatureCollection("projects/bold-listener-471612-n5/assets/zhumadian");
var ROI = zhumadian;
Map.centerObject(ROI, 10);

var startYear = 2023;
var startMonth = 12;
var months = 6;

// MODIS 参考投影
var modisProjection = ee.ImageCollection('MODIS/061/MOD13Q1').filterBounds(ROI).first().projection();

// 1. 去云函数
function maskL8sr(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0)).and(qa.bitwiseAnd(1 << 5).eq(0));
  return image.updateMask(mask);
}

// 2. 指数补全函数 (17个指数，不允许省略)
var addAllIndices = function(img) {
  var b = img.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']).multiply(0.0000275).add(-0.2);
  var Blue = b.select('SR_B2'), Green = b.select('SR_B3'), Red = b.select('SR_B4'), NIR = b.select('SR_B5'), SWIR1 = b.select('SR_B6');

  var ndvi = NIR.subtract(Red).divide(NIR.add(Red)).rename('NDVI');
  var evi = b.expression('2.5 * ((N - R) / (N + 6 * R - 7.5 * B + 1))', {'N':NIR,'R':Red,'B':Blue}).rename('EVI');
  var savi = NIR.subtract(Red).multiply(1.5).divide(NIR.add(Red).add(0.5)).rename('SAVI');
  var osavi = NIR.subtract(Red).divide(NIR.add(Red).add(0.16)).rename('OSAVI');
  var msavi = b.expression('(2*N+1 - sqrt(pow(2*N+1, 2)-8*(N-R)))/2', {'N':NIR,'R':Red}).rename('MSAVI');
  var gndvi = NIR.subtract(Green).divide(NIR.add(Green)).rename('GNDVI');
  var dvi = NIR.subtract(Red).rename('DVI');
  var rvi = NIR.divide(Red).rename('RVI');
  var ci_green = NIR.divide(Green).subtract(1).rename('CI_green');
  var ndmi = NIR.subtract(SWIR1).divide(NIR.add(SWIR1)).rename('NDMI');
  var lswi = NIR.subtract(SWIR1).divide(NIR.add(SWIR1)).rename('LSWI');
  var gvmi = b.expression('((N+0.1)-(S+0.02))/((N+0.1)+(S+0.02))', {'N':NIR,'S':SWIR1}).rename('GVMI');
  var ndwi = Green.subtract(NIR).divide(Green.add(NIR)).rename('NDWI');
  var arvi = b.expression('(N-(2*R-B))/(N+(2*R-B))', {'N':NIR,'R':Red,'B':Blue}).rename('ARVI');
  var vari = Green.subtract(Red).divide(Green.add(Red).subtract(Blue)).rename('VARI');
  var tgi = Green.subtract(Red.multiply(0.39)).subtract(Blue.multiply(0.61)).rename('TGI');
  var exg = Green.multiply(2).subtract(Red).subtract(Blue).rename('ExG');

  return ee.Image([ndvi, evi, savi, osavi, msavi, gndvi, dvi, rvi, ci_green, ndmi, lswi, gvmi, ndwi, arvi, vari, tgi, exg]).multiply(10000).short();
};

// 3. 逐月处理
for (var i = 0; i < months; i++) {
  var date = ee.Date.fromYMD(startYear, startMonth, 1).advance(i, 'month');
  var monthStr = date.format('YYYYMM').getInfo();

  var l8_col = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").filterBounds(ROI).filterDate(date, date.advance(1, 'month')).map(maskL8sr);

  if (l8_col.size().getInfo() === 0) { print(monthStr + ' 无 Landsat 影像，跳过。'); continue; }
  var l8_median = l8_col.median().clip(ROI);
  if (l8_median.bandNames().size().getInfo() === 0) { print(monthStr + ' 去云后无有效像素，跳过。'); continue; }

  var l8_indices = addAllIndices(l8_median);
  var ndvi_raw = l8_indices.select('NDVI').divide(10000);

  // --- 更换为 Gap-Filled GPP 产品 (MOD17A2HGF) ---
  var mod_gpp_col = ee.ImageCollection('MODIS/061/MOD17A2HGF').filterBounds(ROI).filterDate(date, date.advance(1, 'month')).select('Gpp');

  var gpp_img;
  if (mod_gpp_col.size().getInfo() > 0) {
      gpp_img = mod_gpp_col.median().multiply(0.0125); // 还原系数
  } else {
      print(monthStr + ' MODIS GPP (GF) 仍无数据。'); continue;
  }

  var mod_lst = ee.ImageCollection('MODIS/061/MOD11A2').filterBounds(ROI).filterDate(date, date.advance(1, 'month')).select('LST_Day_1km').median().multiply(0.02);

  // 融合逻辑
  var l8_aggr = ndvi_raw.reproject({crs: modisProjection, scale: 250}).reduceResolution(ee.Reducer.mean(), false, 1024);
  var ratio = gpp_img.divide(l8_aggr.add(0.001)).clamp(0.1, 10);
  var ratio_30m = ratio.resample('bilinear').reproject({crs: 'EPSG:4326', scale: 30});

  var final_GPP = ndvi_raw.multiply(ratio_30m).multiply(10000).short().rename('GPP');
  var final_LST = mod_lst.resample('bilinear').reproject({crs: 'EPSG:4326', scale: 30}).multiply(100).short().rename('LST');

  var exportImage = l8_indices.addBands([final_GPP, final_LST]);

  // 导出
  Export.image.toDrive({
    image: exportImage,
    description: 'zhumadian_' + monthStr,
    folder: 'zhumadian_L8_Final',
    region: ROI.geometry(),
    scale: 30,
    crs: 'EPSG:32649',
    maxPixels: 1e13
  });

// ===================== 优化后的可视化对比 =====================
  if (monthStr === '202312' || monthStr === '202404') {

    // 1. 修正后的真彩色：调大 max 到 3000-4000，增加 gamma 解决大面积发白
    Map.addLayer(l8_median, {
      bands: ['SR_B4', 'SR_B3', 'SR_B2'],
      min: 0,
      max: 3500, // 调大这个值解决全白问题
      gamma: 1.4
    }, monthStr + ' 1-真彩色(看底图)');

    // 2. 强力推荐：近红外假彩色 (B5, B4, B3) —— 找回“小格子”的关键
    // 在这个模式下，小麦是红色的，田埂和道路是青色的，对比度最高
    Map.addLayer(l8_median, {
      bands: ['SR_B5', 'SR_B4', 'SR_B3'],
      min: 0,
      max: 5000
    }, monthStr + ' 2-假彩色(找边界/格子)');

    // 3. NDVI 预览
    Map.addLayer(ndvi_raw, {
      min: 0,
      max: 0.9, // 4月调大上限
      palette: ['white', '#99B718', '#056201']
    }, monthStr + ' 3-NDVI');

    // 4. 修正后的 GPP：调大 max 解决过度饱和红色
    Map.addLayer(final_GPP.divide(10000), {
      min: 0,
      max: 1.2, // 4月冬小麦旺盛期 GPP 较高，上限设为 1.2
      palette: ['#FFFFFF', 'blue', 'cyan', 'green', 'yellow', 'red']
    }, monthStr + ' 4-GPP(修正拉伸)');
  }

  print(monthStr + ' 任务已提交。');
}
