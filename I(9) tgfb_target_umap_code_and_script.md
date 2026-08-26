# TGFB 표적 유전자 UMAP 시각화 발표 자료

## 코드

```r
library(ggplot2)
library(patchwork)

# 1. UMAP 차원 축소 좌표 및 주요 표적 유전자 4종의 발현량 데이터 추출
# (Seurat 객체에서 데이터 프레임 형태로 변환하여 시각화 준비)
umap_coords <- as.data.frame(Embeddings(endo_cells, reduction = "umap"))
genes <- c("SERPINE1", "SNAI1", "SNAI2", "TWIST1")
expr_matrix <- LayerData(endo_cells, assay = "RNA", layer = "data")
expr_data <- as.data.frame(t(as.matrix(expr_matrix[genes, ])))

# 2. 시각화를 위해 UMAP 좌표와 유전자 발현량 데이터를 하나의 데이터프레임으로 병합
plot_data <- cbind(umap_coords, expr_data)

# 3. ggplot2를 이용한 각 표적 유전자(TGFB 타겟)의 발현 수준 UMAP 시각화
p1 <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = SERPINE1)) +
  geom_point(size = 0.5) + scale_color_gradient(low = "lightgrey", high = "purple") +
  ggtitle("SERPINE1 (TGFB 타겟)") + theme_minimal()

p2 <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = SNAI1)) +
  geom_point(size = 0.5) + scale_color_gradient(low = "lightgrey", high = "purple") +
  ggtitle("SNAI1 (TGFB 타겟)") + theme_minimal()

p3 <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = SNAI2)) +
  geom_point(size = 0.5) + scale_color_gradient(low = "lightgrey", high = "purple") +
  ggtitle("SNAI2 (TGFB 타겟)") + theme_minimal()

p4 <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = TWIST1)) +
  geom_point(size = 0.5) + scale_color_gradient(low = "lightgrey", high = "purple") +
  ggtitle("TWIST1 (TGFB 타겟)") + theme_minimal()

# 4. patchwork 패키지를 활용한 2x2 형태의 다중 플롯(Multi-plot) 배열 및 출력
(p1 | p2) / (p3 | p4)
```

## 발표 대본

이 코드는 혈관 내피세포 하위 군집 객체인 `endo_cells`를 바탕으로, TGFB 신호전달에 의해 유도될 수 있는 대표적인 표적 유전자들의 발현 분포를 UMAP 위에 시각화하기 위한 코드입니다. 목적은 `SERPINE1`, `SNAI1`, `SNAI2`, `TWIST1`과 같은 TGFB 관련 표적 유전자들이 어떤 세포 집단에서 높게 발현되는지를 확인함으로써, 특정 하위 군집이 TGFB 신호의 영향을 더 강하게 받고 있는지, 그리고 EndoMT와 관련된 전사 프로그램이 어느 위치에서 활성화되는지를 해석하는 것입니다. 먼저 `library(ggplot2)`와 `library(patchwork)`를 불러옵니다. `ggplot2`는 UMAP 기반 산점도 시각화를 만드는 데 사용되고, `patchwork`는 여러 개의 그래프를 하나의 화면에 배열하는 데 사용됩니다.

그다음 `Embeddings(endo_cells, reduction = "umap")`를 이용해 `endo_cells` 객체에 저장된 UMAP 좌표를 추출하고, 이를 데이터프레임으로 변환합니다. 이 좌표는 각 세포가 2차원 UMAP 공간에서 어디에 위치하는지를 나타냅니다. 이후 `genes <- c("SERPINE1", "SNAI1", "SNAI2", "TWIST1")`를 통해 시각화할 네 개의 핵심 유전자를 지정합니다. 이 유전자들은 모두 TGFB 신호전달 또는 EndoMT 관련 해석에서 중요한 표적 유전자로 활용됩니다. 이어서 `LayerData(endo_cells, assay = "RNA", layer = "data")`를 사용하여 RNA assay의 정규화된 발현량 데이터를 가져오고, 그중 지정한 네 개 유전자에 대한 발현량만 추출한 뒤 전치하여 세포별 데이터 형태로 정리합니다. 이후 `cbind()`를 사용하여 UMAP 좌표와 유전자 발현량 정보를 하나의 데이터프레임으로 병합합니다. 이렇게 하면 각 세포에 대해 위치 정보와 표적 유전자 발현 정보가 동시에 포함된 시각화용 데이터가 완성됩니다.

이후 `ggplot()`을 이용해 네 개 유전자 각각에 대한 UMAP 발현 플롯을 생성합니다. `p1`은 `SERPINE1`, `p2`는 `SNAI1`, `p3`는 `SNAI2`, `p4`는 `TWIST1`의 발현량을 나타냅니다. 네 그래프 모두 x축과 y축은 각각 `umap_1`, `umap_2`이고, 각 점은 하나의 세포를 의미하며, 점의 색은 해당 유전자의 발현량을 의미합니다. 발현이 낮은 세포는 연한 회색으로, 발현이 높은 세포는 보라색으로 표시됩니다. 따라서 특정 영역 또는 특정 하위 군집에서 색이 진하게 나타난다면, 그 부분에서 해당 TGFB 표적 유전자가 더 높게 발현되고 있음을 뜻합니다. `SERPINE1`, `SNAI1`, `SNAI2`, `TWIST1`은 TGFB 신호와 관련된 전사 반응이나 EndoMT 진행 과정에서 중요하게 해석될 수 있으므로, 이 시각화는 단순한 발현 확인을 넘어서 세포 상태 변화와 신호 활성의 위치적 분포를 보여 주는 데 의미가 있습니다.

마지막으로 `(p1 | p2) / (p3 | p4)` 구문을 사용하여 네 개의 그래프를 2행 2열 형태로 한 화면에 배치합니다. 이렇게 하면 네 개 유전자의 발현 양상을 동시에 비교할 수 있어, 특정 군집에서 TGFB 관련 프로그램이 일관되게 활성화되는지 더 직관적으로 확인할 수 있습니다. 정리하면, 이 코드는 혈관 내피세포 하위 군집의 UMAP 위에 TGFB 신호전달 관련 표적 유전자 네 개의 발현량을 시각화하여, 특정 세포 집단에서 TGFB 반응성과 EndoMT 관련 전사 프로그램이 강화되어 있는지를 해석하기 위한 코드입니다. 핵심 의의는 하나의 유전자만 보는 것이 아니라, TGFB와 연관된 여러 표적 유전자를 동시에 비교함으로써 세포 상태 변화의 분자적 근거를 보다 설득력 있게 제시할 수 있다는 점에 있습니다.
