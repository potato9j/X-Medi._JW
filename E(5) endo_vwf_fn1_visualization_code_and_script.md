# 혈관 내피세포 하위 군집 마커 시각화 발표 자료

## 코드

```r
library(ggplot2)
library(patchwork)

# 1. UMAP 차원 축소 좌표 및 특정 유전자(VWF, FN1)의 발현량 데이터 추출
# (Seurat 내장 시각화 함수의 버전 충돌 오류를 우회하기 위한 수동 추출)
umap_coords <- as.data.frame(Embeddings(endo_cells, reduction = "umap"))
genes <- c("VWF", "FN1")
expr_matrix <- LayerData(endo_cells, assay = "RNA", layer = "data")
expr_data <- as.data.frame(t(as.matrix(expr_matrix[genes, ])))

# 2. 시각화를 위해 UMAP 좌표와 유전자 발현량 데이터를 하나의 데이터프레임으로 병합
plot_data <- cbind(umap_coords, expr_data)

# 3. ggplot2를 이용한 특정 유전자 발현량 UMAP 시각화
# VWF: 정상 혈관 내피세포(Endothelial cell) 발현 마커
p1 <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = VWF)) +
  geom_point(size = 0.5) + scale_color_gradient(low = "lightgrey", high = "purple") +
  ggtitle("VWF (정상 내피세포 마커)") + theme_minimal()

# FN1: 내피-중간엽 이행(EndoMT) 및 변형된 세포 발현 마커
p2 <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = FN1)) +
  geom_point(size = 0.5) + scale_color_gradient(low = "lightgrey", high = "red") +
  ggtitle("FN1 (EndoMT 마커)") + theme_minimal()

# 4. 두 개의 UMAP 플롯을 나란히 배치하여 비교 시각화
p1 | p2
```

## 발표 대본

이 코드는 앞서 하위 군집 분석을 수행한 `endo_cells` 객체를 바탕으로, 혈관 내피세포 집단 내부에서 `VWF`와 `FN1` 유전자의 발현 분포를 UMAP 상에 시각화하여 세포 상태의 차이를 해석하기 위한 코드입니다. 목적은 정상적인 혈관 내피세포의 특성을 나타내는 마커와, 내피-중간엽 이행 즉 EndoMT와 관련된 변형된 상태를 나타내는 마커를 비교함으로써, 하위 군집 내부에 서로 다른 생물학적 상태가 존재하는지를 확인하는 것입니다. 먼저 `library(ggplot2)`와 `library(patchwork)`를 불러옵니다. `ggplot2`는 UMAP 산점도 기반의 발현 시각화를 만드는 데 사용되고, `patchwork`는 여러 개의 그래프를 한 화면에 나란히 배치하는 데 사용됩니다.

그다음 `Embeddings(endo_cells, reduction = "umap")`를 이용하여 `endo_cells` 객체에 저장된 UMAP 좌표를 추출하고, 이를 데이터프레임으로 변환합니다. 이 좌표는 각 세포가 2차원 UMAP 공간에서 어디에 위치하는지를 나타냅니다. 이후 `genes <- c("VWF", "FN1")`를 통해 시각화할 두 개의 핵심 유전자를 지정합니다. 여기서 `VWF`는 정상적인 혈관 내피세포의 대표적인 마커이고, `FN1`은 내피-중간엽 이행인 EndoMT나 보다 변형된 세포 상태와 관련된 마커입니다. 이어서 `LayerData(endo_cells, assay = "RNA", layer = "data")`를 사용해 RNA assay의 정규화된 발현값을 불러오고, 선택한 두 유전자에 해당하는 발현량만 추출한 뒤 전치하여 세포별 발현 데이터 형태로 변환합니다. 그 다음 `cbind()`를 이용해 UMAP 좌표와 유전자 발현량 데이터를 하나의 데이터프레임으로 병합합니다. 이 과정이 끝나면 각 세포에 대해 위치 정보와 두 유전자의 발현량이 동시에 정리된 시각화용 데이터가 완성됩니다.

이후 `ggplot()`을 사용하여 두 개의 UMAP 발현 플롯을 각각 생성합니다. 첫 번째 그래프인 `p1`은 `VWF` 발현량을 표시하는 그래프로, 정상 혈관 내피세포의 특성이 어느 위치의 세포에서 강하게 나타나는지를 보여 줍니다. 두 번째 그래프인 `p2`는 `FN1` 발현량을 표시하는 그래프로, EndoMT와 연관된 변형된 세포 상태가 UMAP 상 어느 부분에 분포하는지를 보여 줍니다. 두 그래프 모두 x축과 y축은 각각 `umap_1`, `umap_2`이고, 각 점은 하나의 세포를 의미하며, 점의 색은 해당 유전자의 발현량을 나타냅니다. 발현이 낮은 세포는 연한 회색으로, 높은 세포는 지정된 색상으로 표시됩니다. `VWF`는 보라색 계열, `FN1`은 빨간색 계열로 표현하여 두 유전자의 생물학적 의미를 시각적으로 구분하기 쉽게 했습니다. 마지막으로 `p1 | p2`를 사용하여 두 개의 그래프를 가로로 나란히 배치함으로써, 정상 내피세포 마커와 EndoMT 마커의 발현 패턴을 직접 비교할 수 있도록 구성합니다. 정리하면, 이 코드는 혈관 내피세포 하위 군집의 UMAP 위에 `VWF`와 `FN1`의 발현량을 시각화하여, 세포들이 정상 내피 특성을 유지하는지 또는 EndoMT와 같은 변화된 상태를 보이는지를 비교 해석하기 위한 코드입니다. 핵심 의의는 동일한 세포 집단 내부에서도 서로 다른 기능적 상태가 존재할 수 있음을 마커 발현 패턴을 통해 보다 명확하게 보여 준다는 점에 있습니다.
