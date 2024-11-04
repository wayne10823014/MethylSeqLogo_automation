// 定義顏色方案
const COLOR_SCHEME = {
    'G': 'orange',
    'A': 'green',
    'C': 'deepskyblue',
    'T': 'red'
};

// 函數：將數據轉換為 LogoJS 的格式
function convertDataToLogoJSFormat(baseHeights) {
    return baseHeights.map((baseHeight) => {
        return baseHeight.map(([base, height]) => ({
            base: base,
            height: height,
            color: COLOR_SCHEME[base]
        }));
    });
}

// 模擬數據：這裡是轉換後的數據（可以根據需要替換為實際數據）
const baseHeights = [
    [['A', 2.5], ['T', 1.5], ['C', 1.0], ['G', 0.5]],
    [['C', 2.0], ['G', 2.0], ['A', 1.0], ['T', 1.0]],
    [['G', 2.5], ['C', 1.5], ['T', 1.0], ['A', 0.5]],
    // 其他位置數據...
];

// 創建 LogoJS 的選項
const options = {
    data: convertDataToLogoJSFormat(baseHeights),
    height: 300, // 圖像高度
    width: 800,  // 圖像寬度
    padding: 10, // 每個標誌之間的間距
    fontSize: 14 // 字體大小
};

// 使用 LogoJS 繪製序列標誌
window.onload = function() {
    renderLogo(document.getElementById('logo-container'), options);
};
