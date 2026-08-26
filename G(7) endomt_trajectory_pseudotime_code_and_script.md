# EndoMT 궤적 및 가상 시간 분석 발표 자료

## 코드

```r
# =====================================================================
# 세포 분화 궤적(Trajectory) 및 가상 시간(Pseudotime) 분석
# =====================================================================

# 1. 필수 패키지 로드 (궤적 분석 및 색상 팔레트 도구)
library(slingshot)
library(RColorBrewer)

# 2. Slingshot 분석을 위한 UMAP 좌표 및 클러스터 메타데이터 추출
umap_coords <- Embeddings(endo_cells, reduction = "umap")
cluster_labels <- Idents(endo_cells)

# 3. [핵심 복원] Slingshot 알고리즘 실행 (sds 객체 생성)
# 정상 내피세포인 9번 클러스터를 시작점(start.clus)으로 지정하여 분화 궤적 계산
sds <- slingshot(data = umap_coords, clusterLabels = cluster_labels, start.clus = "9")

# 4. 계산된 궤적에서 가상 시간(Pseudotime) 값 추출
pt <- slingPseudotime(sds)[, 1]

# 5. 가상 시간의 흐름을 나타낼 연속적인 색상 팔레트 설정 (Spectral)
colors <- colorRampPalette(brewer.pal(11, 'Spectral')[-6])(100)
plotcol <- colors[cut(pt, breaks = 100)]
plotcol[is.na(plotcol)] <- "lightgrey" # 궤적에 포함되지 않은 세포는 회색 처리

# 6. 최종 궤적 시각화 (UMAP 상에 Slingshot 곡선 투영)
plot(umap_coords, col = plotcol, pch = 16, cex = 0.5,
     main = "EndoMT Trajectory (정상 9번 -> EndoMT 12/13번)",
     xlab = "UMAP_1", ylab = "UMAP_2")
lines(SlingshotDataSet(sds), lwd = 3, col = 'black')
```

## 발표 대본

이 코드는 혈관 내피세포 하위 군집 객체인 `endo_cells`를 대상으로, 세포 상태 변화의 방향성을 확인하기 위해 궤적 분석과 가상 시간 분석을 수행하는 코드입니다. 목적은 정상 내피세포로 해석되는 클러스터 9번을 출발점으로 두고, EndoMT 성격을 보이는 클러스터 12번과 13번 방향으로 세포 상태가 어떻게 이어지는지를 연속적인 흐름으로 시각화하는 것입니다. 즉, 이 분석은 단순히 서로 다른 군집이 존재한다는 사실을 보여주는 데서 그치지 않고, 세포들이 하나의 상태에서 다른 상태로 변화하는 과정이 추정 가능한지를 확인하는 데 의미가 있습니다. 먼저 `library(slingshot)`과 `library(RColorBrewer)`를 불러옵니다. `slingshot`은 단일세포 데이터에서 분화 궤적과 가상 시간을 계산하는 데 사용되는 대표적인 패키지이고, `RColorBrewer`는 가상 시간의 흐름을 색으로 표현하기 위한 팔레트를 제공하는 도구입니다.

그다음 `Embeddings(endo_cells, reduction = "umap")`를 사용하여 `endo_cells` 객체에 저장된 UMAP 좌표를 추출하고, `Idents(endo_cells)`를 통해 현재 세포들의 군집 정보를 가져옵니다. 여기서 UMAP 좌표는 세포들 사이의 저차원 공간상 위치 정보를 의미하고, 군집 라벨은 Slingshot이 궤적을 연결할 때 참고하는 집단 정보입니다. 이후 핵심 단계인 `sds <- slingshot(data = umap_coords, clusterLabels = cluster_labels, start.clus = "9")`를 실행합니다. 이 코드는 Slingshot 알고리즘을 이용하여 세포 상태 전이 경로를 계산하는 부분이며, `start.clus = "9"`는 9번 클러스터를 시작점으로 지정한다는 뜻입니다. 즉, 분석자는 정상 내피세포로 해석되는 9번 클러스터를 출발 상태로 가정하고, 그로부터 다른 세포 상태들로 이어지는 가능한 생물학적 경로를 추정하도록 설정한 것입니다.

이후 `pt <- slingPseudotime(sds)[, 1]`을 사용하여 계산된 궤적에서 각 세포의 가상 시간 값을 추출합니다. 가상 시간은 실제 시간 측정값이 아니라, 세포가 추정된 분화 또는 상태 변화 경로상에서 어느 정도 앞이나 뒤에 위치하는지를 상대적으로 나타내는 값입니다. 따라서 값이 작은 세포는 시작점에 더 가까운 상태이고, 값이 큰 세포는 더 진행된 상태로 해석할 수 있습니다. 그다음 `colorRampPalette(brewer.pal(11, 'Spectral')[-6])(100)`을 사용하여 가상 시간의 연속적인 변화를 표현할 색상 팔레트를 생성합니다. 이후 `plotcol <- colors[cut(pt, breaks = 100)]`를 통해 각 세포의 가상 시간 값에 맞는 색을 배정하고, `plotcol[is.na(plotcol)] <- "lightgrey"`를 사용해 궤적 계산에 포함되지 않은 세포는 회색으로 표시합니다. 이 설정을 통해 궤적에 포함된 세포와 그렇지 않은 세포를 시각적으로 구분할 수 있습니다.

마지막으로 `plot()`을 이용해 UMAP 좌표 위에 각 세포를 점으로 표시하고, 점의 색을 가상 시간에 따라 다르게 나타냅니다. 여기서 `main = "EndoMT Trajectory (정상 9번 -> EndoMT 12/13번)"`는 그래프 제목이고, `xlab`, `ylab`은 UMAP 축 이름을 표시합니다. 이어서 `lines(SlingshotDataSet(sds), lwd = 3, col = 'black')`를 사용하여 Slingshot이 계산한 궤적 곡선을 UMAP 위에 검은색 선으로 덧그립니다. 이 곡선은 정상 내피세포 상태에서 EndoMT 성격의 세포 상태로 이어지는 추정 경로를 의미합니다. 정리하면, 이 코드는 혈관 내피세포 하위 군집에서 정상 상태로 해석되는 클러스터를 시작점으로 설정한 뒤, Slingshot 알고리즘을 이용해 세포 상태 변화의 연속적인 경로와 각 세포의 가상 시간을 계산하고, 이를 UMAP 위에 시각화하는 분석 코드입니다. 핵심 의의는 내피세포와 EndoMT 세포가 단순히 떨어진 별개의 집단이 아니라, 하나의 연속적인 상태 변화 과정으로 연결될 가능성을 시각적으로 제시할 수 있다는 점에 있습니다.
