# UMAP 마커 유전자 시각화 발표 자료

## 코드

```r
library(ggplot2)
library(patchwork)

# 1. 시각화를 위한 UMAP 좌표 데이터 추출 (오류 우회를 위한 수동 추출)
umap_coords <- as.data.frame(Embeddings(merged_seurat, reduction = "umap"))

# 2. 대식세포(Macrophage) 및 혈관 내피세포(Endothelial) 핵심 마커 유전자 발현량 추출
genes <- c("CD163", "CD68", "PECAM1", "CLDN5")
expr_matrix <- LayerData(merged_seurat, assay = "RNA", layer = "data")
expr_data <- as.data.frame(t(as.matrix(expr_matrix[genes, ])))

# 3. UMAP 좌표와 유전자 발현량 데이터를 하나의 데이터 프레임으로 병합
plot_data <- cbind(umap_coords, expr_data)

# 4. 각 세포 군집별 마커 유전자 발현량 시각화 (ggplot2 활용)
p1 <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = CD163)) +
  geom_point(size = 0.1) + scale_color_gradient(low = "lightgrey", high = "purple") +
  ggtitle("CD163 (Macrophage 마커)") + theme_minimal()

p2 <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = CD68)) +
  geom_point(size = 0.1) + scale_color_gradient(low = "lightgrey", high = "purple") +
  ggtitle("CD68 (Macrophage 마커)") + theme_minimal()

p3 <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = PECAM1)) +
  geom_point(size = 0.1) + scale_color_gradient(low = "lightgrey", high = "red") +
  ggtitle("PECAM1 (Endothelial 마커)") + theme_minimal()

p4 <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = CLDN5)) +
  geom_point(size = 0.1) + scale_color_gradient(low = "lightgrey", high = "red") +
  ggtitle("CLDN5 (Endothelial 마커)") + theme_minimal()

# 5. 4개의 그래프를 2x2 배열로 한 화면에 모아서 출력
(p1 | p2) / (p3 | p4)
```

## 발표 대본

이 코드는 앞서 생성된 UMAP 결과를 바탕으로, 특정 마커 유전자의 발현량을 시각화하여 각 세포 군집의 생물학적 정체성을 해석하기 위한 코드입니다. 구체적으로는 대식세포를 대표하는 마커인 `CD163`, `CD68`과 혈관 내피세포를 대표하는 마커인 `PECAM1`, `CLDN5`의 발현 분포를 UMAP 위에 표시하여, 어떤 군집이 대식세포 성격을 가지는지, 어떤 군집이 혈관 내피세포 성격을 가지는지를 시각적으로 확인하는 것이 목적입니다. 먼저 `library(ggplot2)`와 `library(patchwork)`를 사용하여 시각화에 필요한 패키지를 불러옵니다. `ggplot2`는 UMAP 산점도를 그리는 데 사용되고, `patchwork`는 여러 개의 그래프를 한 화면에 배열하는 데 사용됩니다.

그다음 `Embeddings(merged_seurat, reduction = "umap")`를 이용해 Seurat 객체에 저장된 UMAP 좌표를 추출하고, 이를 데이터프레임으로 변환합니다. 이 단계는 각 세포가 UMAP 공간상 어디에 위치하는지를 가져오는 과정입니다. 이후 `genes <- c("CD163", "CD68", "PECAM1", "CLDN5")`로 시각화할 핵심 마커 유전자를 지정하고, `LayerData(merged_seurat, assay = "RNA", layer = "data")`를 사용하여 RNA assay의 정규화된 발현값을 불러옵니다. 그다음 선택한 네 개 유전자에 대한 발현량만 추출한 뒤, 전치를 통해 세포별 발현 정보 형태로 바꾸고 데이터프레임으로 정리합니다. 이어서 `cbind()`를 사용하여 UMAP 좌표와 유전자 발현량 데이터를 하나의 데이터프레임으로 병합합니다. 이렇게 하면 각 세포에 대해 위치 정보와 마커 발현 정보가 동시에 정리되어 이후 시각화에 바로 사용할 수 있는 구조가 완성됩니다.

이후 `ggplot()`을 이용해 네 개의 마커 유전자에 대한 시각화를 각각 생성합니다. `p1`은 `CD163`, `p2`는 `CD68`, `p3`는 `PECAM1`, `p4`는 `CLDN5`의 발현량을 UMAP 위에 표시한 그래프입니다. 이때 x축과 y축은 각각 `umap_1`, `umap_2`이고, 점의 색은 해당 유전자의 발현량을 의미합니다. 발현이 낮은 세포는 연한 회색으로, 높은 세포는 지정된 색으로 나타납니다. `CD163`와 `CD68`은 대식세포 마커이기 때문에 보라색 계열로 표현하였고, `PECAM1`과 `CLDN5`는 혈관 내피세포 마커이기 때문에 빨간색 계열로 표현하였습니다. 따라서 특정 군집에서 `CD163`와 `CD68`이 동시에 높게 보이면 대식세포 성격을 가진 군집으로 해석할 수 있고, `PECAM1`과 `CLDN5`가 높게 보이면 혈관 내피세포 성격을 가진 군집으로 해석할 수 있습니다.

마지막으로 `(p1 | p2) / (p3 | p4)` 구문을 사용하여 네 개의 그래프를 2행 2열 배열로 한 화면에 배치합니다. 이렇게 하면 대식세포 마커 두 개와 혈관 내피세포 마커 두 개를 동시에 비교할 수 있어서, 각 군집의 세포 유형을 더 직관적으로 해석할 수 있습니다. 정리하면, 이 코드는 이미 계산된 UMAP 좌표를 기반으로 주요 마커 유전자 발현량을 시각화하여, 세포 군집이 어떤 세포 유형에 해당하는지를 해석하는 과정입니다. 핵심 의의는 단순히 군집의 위치만 보는 것이 아니라, 마커 유전자 발현 패턴을 함께 확인함으로써 각 군집의 생물학적 특성을 보다 명확하게 설명할 수 있다는 점에 있습니다.
