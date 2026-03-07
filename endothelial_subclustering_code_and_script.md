# 혈관 내피세포 하위 군집 분석 발표 자료

## PPT에 게시할 코드

```r
# =====================================================================
# 혈관 내피세포(Endothelial cells) 하위 군집 분석 (Sub-clustering)
# =====================================================================

# 1. [핵심 추가!] 전체 데이터에서 9, 12, 13, 18번 클러스터(혈관세포)만 추출하여 새로운 객체 생성
endo_cells <- subset(merged_seurat, idents = c(9, 12, 13, 18))

# 2. 하위 분석의 기준을 '하모니 통합본'이 아닌 '원본 발현량(RNA)'으로 초기화
DefaultAssay(endo_cells) <- "RNA"

# 3. Seurat V5 에러 방지를 위해 데이터 형식을 최신(Assay5)으로 변환
endo_cells[["RNA"]] <- as(object = endo_cells[["RNA"]], Class = "Assay5")

# 4. 추출된 혈관세포들만을 대상으로 다시 정규화 및 PCA 진행
endo_cells <- NormalizeData(endo_cells)
endo_cells <- FindVariableFeatures(endo_cells)
endo_cells <- ScaleData(endo_cells)
endo_cells <- RunPCA(endo_cells, verbose = FALSE)

# 5. 혈관세포들만의 새로운 2차원 UMAP 생성 (형태 고정을 위해 seed 42 사용)
endo_cells <- RunUMAP(endo_cells, dims = 1:20, seed.use = 42)

# 6. 최종 추출된 혈관세포 UMAP 시각화 (깔끔한 클래식 테마 적용)
library(ggplot2)
DimPlot(endo_cells, reduction = "umap", label = TRUE, pt.size = 0.5) + theme_classic()
```

## 발표 대본

이 코드는 전체 단일세포 RNA-seq 데이터에서 혈관 내피세포로 해석된 특정 클러스터만 따로 추출한 뒤, 그 세포들만 대상으로 다시 차원 축소와 시각화를 수행하여 혈관 내피세포 내부의 세부적인 하위 구조를 분석하기 위한 코드입니다. 즉, 전체 데이터에서 이미 혈관세포로 추정된 집단을 다시 분리해서, 그 안에 서로 다른 아형이나 상태 차이가 존재하는지를 확인하는 것이 목적입니다. 먼저 `subset(merged_seurat, idents = c(9, 12, 13, 18))`를 사용하여 전체 데이터에서 9번, 12번, 13번, 18번 클러스터에 해당하는 세포만 추출하고, 이를 `endo_cells`라는 새로운 Seurat 객체로 저장합니다. 이 단계는 전체 세포 집단이 아니라 혈관 내피세포 집단만 별도로 떼어내어 보다 정밀한 하위 분석을 수행하기 위한 출발점입니다.

그다음 `DefaultAssay(endo_cells) <- "RNA"`를 통해 하위 분석의 기준을 Harmony로 통합된 결과가 아니라 원본 RNA 발현값으로 다시 맞춥니다. 이는 하위 군집 분석에서는 통합본보다 원래의 발현 신호를 기준으로 세포 간 차이를 다시 계산하겠다는 의미입니다. 이어서 `endo_cells[["RNA"]] <- as(object = endo_cells[["RNA"]], Class = "Assay5")`를 사용하여 RNA assay를 Seurat V5의 최신 형식인 `Assay5`로 변환합니다. 이 과정은 분석 자체의 의미를 바꾸는 단계라기보다, Seurat V5 환경에서 발생할 수 있는 데이터 구조 관련 오류를 방지하기 위한 호환성 조정 단계입니다.

이후에는 추출된 혈관 내피세포 집단만을 대상으로 다시 전처리를 수행합니다. `NormalizeData()`는 세포별 총 발현량 차이를 보정하여 비교 가능한 상태로 만들고, `FindVariableFeatures()`는 이 하위 집단 내부에서 세포 간 차이를 잘 설명하는 고변동 유전자를 다시 선택합니다. `ScaleData()`는 유전자별 발현값의 스케일을 맞추어 특정 유전자가 과도하게 분석을 지배하지 않도록 조정하며, `RunPCA()`는 혈관 내피세포들만의 발현 패턴을 바탕으로 저차원 공간을 새롭게 계산합니다. 이 과정은 전체 데이터에서 한 번 수행했던 전처리를 그대로 재사용하는 것이 아니라, 혈관세포 집단 안에서만 다시 최적화된 구조를 계산한다는 점에서 중요합니다.

그다음 `RunUMAP(endo_cells, dims = 1:20, seed.use = 42)`를 사용하여 혈관 내피세포들만의 새로운 UMAP을 생성합니다. 여기서 `dims = 1:20`은 PCA의 1번부터 20번 차원까지를 사용한다는 의미이고, `seed.use = 42`는 UMAP 결과가 실행할 때마다 크게 달라지지 않도록 형태를 고정하기 위한 난수 시드 설정입니다. 즉, 혈관세포 집단 내부의 구조를 2차원 공간에서 안정적으로 시각화하기 위한 단계입니다. 마지막으로 `library(ggplot2)`를 불러오고, `DimPlot(endo_cells, reduction = "umap", label = TRUE, pt.size = 0.5) + theme_classic()`를 실행하여 최종 UMAP을 시각화합니다. 이때 `label = TRUE`는 군집 레이블을 표시하는 설정이고, `pt.size = 0.5`는 점 크기를 조정하는 옵션이며, `theme_classic()`은 불필요한 배경 요소를 줄인 클래식 스타일을 적용하여 그래프를 더 깔끔하게 보이도록 합니다. 정리하면, 이 코드는 전체 데이터에서 혈관 내피세포로 추정된 클러스터만 따로 추출한 뒤, 원본 RNA 발현값을 기준으로 다시 정규화와 PCA를 수행하고, 새로운 UMAP을 생성하여 혈관 내피세포 내부의 세부적인 하위 군집 구조를 탐색하기 위한 분석 코드입니다. 핵심 의의는 전체 데이터 수준의 큰 군집 구조를 넘어서, 특정 세포 계통 내부의 이질성과 세부 아형을 보다 정밀하게 확인할 수 있다는 점에 있습니다.
