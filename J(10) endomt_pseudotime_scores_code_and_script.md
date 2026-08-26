# EndoMT 가상 시간 기반 점수 변화 시각화 발표 자료

## 코드

```r
library(ggplot2)
library(patchwork)
library(slingshot)

# 1. 분화 궤적 분석(Slingshot) 재실행 및 가상 시간(Pseudotime) 데이터 산출
umap_coords <- Embeddings(endo_cells, "umap")
cell_clusters <- Idents(endo_cells)
sds <- slingshot(data = umap_coords, clusterLabels = cell_clusters, start.clus = "9")
pt <- slingPseudotime(sds)[, 1]

# 2. 정규화된 유전자 발현량 데이터 추출 및 세포 상태 지표(Score) 산출
expr_matrix <- LayerData(endo_cells, assay = "RNA", layer = "data")

# [Endothelial Score]: 내피세포 마커 3종(PECAM1, VWF, CDH5)의 평균 발현 수준 계산
endo_genes <- c("PECAM1", "VWF", "CDH5")
endo_score <- colMeans(as.matrix(expr_matrix[endo_genes, ]))

# [Mesenchymal Score]: 중간엽 세포 마커 4종(FN1, ACTA2, TAGLN, COL1A1)의 평균 발현 수준 계산
mes_genes <- c("FN1", "ACTA2", "TAGLN", "COL1A1")
mes_score <- colMeans(as.matrix(expr_matrix[mes_genes, ]))

# 3. 시각화 및 통계 분석을 위한 통합 데이터 프레임(plot_df) 구축
plot_df <- data.frame(
  Pseudotime = pt,
  Endo_Score = endo_score,
  Mes_Score = mes_score,
  VWF = as.numeric(expr_matrix["VWF", ]),
  FN1 = as.numeric(expr_matrix["FN1", ])
)

# 분석 유효 범위(Trajectory 내 존재 세포)를 확보하기 위한 결측치(NA) 정제
plot_df <- na.omit(plot_df)

# 4. 가상 시간의 흐름에 따른 형질 전환 양상 다중 패널 시각화
# 가독성 증대를 위해 LOESS(Locally Estimated Scatterplot Smoothing) 회귀 곡선 적용
p_score_endo <- ggplot(plot_df, aes(x = Pseudotime, y = Endo_Score)) + 
  geom_smooth(color = "blue", size = 1.5) + theme_minimal() + ggtitle("Endothelial Score ( ↓ )")

p_score_mes <- ggplot(plot_df, aes(x = Pseudotime, y = Mes_Score)) + 
  geom_smooth(color = "red", size = 1.5) + theme_minimal() + ggtitle("Mesenchymal Score ( ↑ )")

p_vwf <- ggplot(plot_df, aes(x = Pseudotime, y = VWF)) + 
  geom_smooth(color = "blue") + theme_minimal() + ggtitle("VWF Expression ( ↓ )")

p_fn1 <- ggplot(plot_df, aes(x = Pseudotime, y = FN1)) + 
  geom_smooth(color = "red") + theme_minimal() + ggtitle("FN1 Expression ( ↑ )")

# 5. 지표별 변화 추이 비교를 위한 2x2 레이아웃 배치 및 출력
(p_score_endo | p_score_mes) / (p_vwf | p_fn1)
```

## 발표 대본

이 코드는 혈관 내피세포 하위 군집 객체인 `endo_cells`를 대상으로, 가상 시간의 흐름에 따라 내피세포 특성과 중간엽 특성이 어떻게 변화하는지를 정량적으로 시각화하기 위한 코드입니다. 목적은 Slingshot을 이용해 계산한 pseudotime을 기준으로, 정상 내피세포 성격을 나타내는 점수는 감소하는지, 중간엽 성격을 나타내는 점수는 증가하는지, 그리고 대표 마커 유전자인 `VWF`와 `FN1`의 발현이 실제로 같은 방향의 변화를 보이는지를 함께 확인하는 것입니다. 먼저 `library(ggplot2)`, `library(patchwork)`, `library(slingshot)`을 불러옵니다. `ggplot2`는 점수 변화 시각화를 위한 그래프를 생성하는 데 사용되고, `patchwork`는 여러 그래프를 한 화면에 배치하는 데 사용되며, `slingshot`은 세포 상태 변화의 궤적과 가상 시간을 계산하는 데 사용됩니다.

그다음 `Embeddings(endo_cells, "umap")`를 이용하여 UMAP 좌표를 추출하고, `Idents(endo_cells)`를 통해 각 세포의 군집 정보를 가져옵니다. 이어서 `slingshot(data = umap_coords, clusterLabels = cell_clusters, start.clus = "9")`를 실행하여 궤적 분석을 다시 수행합니다. 여기서 `start.clus = "9"`는 정상 내피세포로 해석되는 9번 클러스터를 시작점으로 설정한 것이고, 이 기준을 바탕으로 다른 세포 상태로 이어지는 연속적인 경로를 추정합니다. 이후 `slingPseudotime(sds)[, 1]`을 사용해 각 세포의 pseudotime 값을 추출합니다. 이 값은 실제 시간이 아니라, 추정된 상태 변화 경로상에서 세포가 어느 정도 진행된 위치에 있는지를 나타내는 상대적 지표입니다.

이후에는 정규화된 RNA 발현값을 불러오고, 이를 바탕으로 세포 상태를 대표하는 점수를 계산합니다. `LayerData(endo_cells, assay = "RNA", layer = "data")`는 RNA assay의 정규화된 발현량을 가져오는 부분입니다. 그다음 `endo_genes <- c("PECAM1", "VWF", "CDH5")`를 정의한 뒤 `colMeans(as.matrix(expr_matrix[endo_genes, ]))`를 사용하여 각 세포의 Endothelial Score를 계산합니다. 이는 정상 혈관 내피세포 특성을 나타내는 대표 마커 세 개의 평균 발현값입니다. 같은 방식으로 `mes_genes <- c("FN1", "ACTA2", "TAGLN", "COL1A1")`를 정의하고, 네 개 중간엽 마커의 평균 발현값을 계산하여 Mesenchymal Score를 얻습니다. 즉, 이 단계는 개별 유전자 하나가 아니라 기능적으로 관련된 마커 묶음을 평균화하여 세포 상태를 더 안정적으로 수치화하는 과정입니다.

그다음 `plot_df <- data.frame(...)`를 사용하여 pseudotime, Endothelial Score, Mesenchymal Score, 그리고 대표 유전자인 `VWF`, `FN1`의 발현값을 하나의 데이터프레임으로 통합합니다. 이후 `plot_df <- na.omit(plot_df)`를 통해 pseudotime이 없는 세포, 즉 궤적 분석에 포함되지 않은 세포를 제거하여 유효한 분석 범위만 남깁니다. 이렇게 준비된 데이터프레임은 가상 시간에 따른 연속적 변화 양상을 시각화하는 입력 자료가 됩니다.

마지막으로 `ggplot()`과 `geom_smooth()`를 이용해 네 개의 변화를 각각 시각화합니다. `p_score_endo`는 pseudotime에 따른 Endothelial Score 변화를 보여 주고, `p_score_mes`는 Mesenchymal Score의 변화를 보여 줍니다. 또 `p_vwf`는 대표 내피 마커인 `VWF` 발현의 변화, `p_fn1`은 대표 중간엽 마커인 `FN1` 발현의 변화를 나타냅니다. 여기서 `geom_smooth()`는 LOESS 기반의 부드러운 회귀 곡선을 그려서 전체적인 추세를 보기 쉽게 만드는 역할을 합니다. 파란색 그래프는 내피 특성과 관련된 감소 방향을, 빨간색 그래프는 중간엽 특성과 관련된 증가 방향을 직관적으로 보여 주도록 구성되어 있습니다. 마지막으로 `(p_score_endo | p_score_mes) / (p_vwf | p_fn1)`를 사용하여 네 개의 그래프를 2행 2열로 배치함으로써, pseudotime에 따른 상태 전환 양상을 한 화면에서 종합적으로 비교할 수 있도록 합니다. 정리하면, 이 코드는 Slingshot으로 계산한 가상 시간을 기준으로 내피세포 특성과 중간엽 특성의 점수, 그리고 대표 유전자 발현의 변화를 동시에 시각화하여, EndoMT 과정이 연속적인 상태 변화로 진행된다는 점을 정량적으로 제시하기 위한 코드입니다. 핵심 의의는 UMAP 기반의 위치 정보만 보는 데서 나아가, pseudotime 축을 따라 세포 상태가 실제로 어떻게 전환되는지를 점수와 마커 발현 수준으로 함께 설명할 수 있다는 점에 있습니다.
