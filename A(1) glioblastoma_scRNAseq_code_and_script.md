# Glioblastoma scRNA-seq 분석 발표 자료

## 코드

```r
library(Seurat)
library(harmony)

# 1. 전처리된 데이터 불러오기 및 품질 관리(QC) 필터링 적용
# (유전자 발현량이 너무 적거나 많지 않고, 미토콘드리아 비율이 15% 미만인 건강한 세포만 추출)
merged_seurat <- readRDS("Glioblastoma_merged_initial.rds")
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 15)

# 2. 데이터 정규화, 고변동 유전자(Variable Features) 추출 및 선형 차원 축소(PCA)
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat, verbose = FALSE)

# ---------------------------------------------------------
# 3. 하모니(Harmony) 알고리즘을 이용한 배치 효과(Batch effect) 보정
# ---------------------------------------------------------
# 3-1. 하모니 입력용으로 PCA 결과값과 메타데이터(세포 정보) 추출
pca_matrix <- Embeddings(merged_seurat, reduction = "pca")
meta_data <- merged_seurat@meta.data

# 3-2. 샘플 출신지(orig.ident)별 기술적 오차를 통합하기 위해 하모니 알고리즘 실행
harmony_matrix <- HarmonyMatrix(
  data_mat = pca_matrix,
  meta_data = meta_data,
  vars_use = "orig.ident",
  do_pca = FALSE
)

# 3-3. 통합 및 보정된 하모니 결과값을 Seurat 객체 내에 새로운 차원으로 저장
merged_seurat[["harmony"]] <- CreateDimReducObject(embeddings = harmony_matrix, key = "harmony_", assay = DefaultAssay(merged_seurat))
# ---------------------------------------------------------

# 4. 그래프 기반 군집화(Clustering) 및 비선형 차원 축소 (UMAP) 적용
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5)

# 5. 최종 군집화된 UMAP 시각화
DimPlot(merged_seurat, reduction = "umap", label = TRUE, pt.size = 0.1)
```

## 발표 대본

안녕하세요. 제가 설명드릴 코드는 **Glioblastoma 단일세포 RNA-seq 데이터**를 분석하는 전체 파이프라인입니다. 이 코드의 목적은 단일세포 데이터를 불러온 뒤, 품질이 낮은 세포를 제거하고, 샘플 간에 발생하는 기술적 차이인 **배치 효과를 보정한 다음**, 유사한 세포끼리 군집화해서 최종적으로 시각화하는 것입니다.

먼저 코드의 첫 부분에서는 `Seurat`와 `harmony` 패키지를 불러옵니다. `Seurat`는 단일세포 데이터 분석에서 가장 많이 사용하는 패키지로, 전처리, 차원 축소, 군집화, 시각화를 담당합니다. `harmony`는 여러 샘플을 함께 분석할 때 생기는 배치 효과를 줄이는 데 사용됩니다.

그다음 `readRDS()`를 사용해서 전처리된 Seurat 객체를 불러옵니다. 이후 `subset()`을 이용해서 **품질 관리, 즉 QC 필터링**을 수행합니다. 여기서 `nFeature_RNA > 500` 조건은 유전자 검출 수가 너무 적은 세포를 제거하기 위한 것이고, `nFeature_RNA < 6000` 조건은 비정상적으로 많은 유전자가 검출된 세포를 제거해서 doublet 가능성을 줄이기 위한 것입니다. 또 `percent.mt < 15` 조건은 미토콘드리아 유전자 비율이 높은 세포를 제거하기 위한 기준입니다. 따라서 이 단계는 분석에 적합한 세포만 남기는 과정이라고 볼 수 있습니다.

이후에는 데이터를 본격적으로 분석 가능한 형태로 정리합니다. `NormalizeData()`는 세포마다 총 발현량이 다른 문제를 보정해서 서로 비교할 수 있도록 만들어 줍니다. `FindVariableFeatures()`는 세포 간 차이를 잘 설명하는 **고변동 유전자**를 선택하는 단계입니다. `ScaleData()`는 유전자 발현값의 스케일을 맞춰서 특정 유전자만 과도하게 영향을 주지 않도록 조정합니다. 그리고 `RunPCA()`는 고차원 유전자 발현 데이터를 저차원 공간으로 줄여서 데이터의 핵심 변동 구조를 요약합니다. 이 PCA 결과는 뒤에 나오는 배치 효과 보정과 군집화의 기반이 됩니다.

이 코드에서 가장 핵심적인 부분은 **Harmony를 이용한 배치 효과 보정**입니다. 먼저 `Embeddings()`를 이용해서 PCA 좌표를 추출하고, `meta.data`에서 세포별 메타데이터를 가져옵니다. 그 다음 `HarmonyMatrix()`를 실행하는데, 여기서 `vars_use = "orig.ident"`로 설정해서 각 세포가 어느 샘플에서 왔는지를 기준으로 배치 효과를 보정합니다. 즉, 샘플이 다르기 때문에 인위적으로 멀어 보이는 세포들을 조정해서, 실제 생물학적으로 유사한 세포들이 더 가깝게 보이도록 만드는 과정입니다. 그리고 이렇게 보정된 결과는 `CreateDimReducObject()`를 통해 Seurat 객체 안에 `harmony`라는 새로운 차원 축소 결과로 저장됩니다.

그다음 단계에서는 이 Harmony 결과를 기반으로 군집화를 진행합니다. `RunUMAP()`은 보정된 데이터를 2차원 공간에 시각화하는 단계이고, `FindNeighbors()`는 세포 간 유사도 그래프를 만드는 단계입니다. 이후 `FindClusters()`는 이 그래프를 바탕으로 비슷한 세포끼리 같은 군집으로 묶습니다. 여기서 `dims = 1:30`은 Harmony 차원 1번부터 30번까지를 사용했다는 의미이고, `resolution = 0.5`는 군집을 어느 정도 세분화할지를 정하는 값입니다.

마지막으로 `DimPlot()`을 이용해서 UMAP 위에 최종 군집 결과를 시각화합니다. 이때 `label = TRUE`로 설정했기 때문에 각 군집 번호가 함께 표시되고, `pt.size = 0.1`로 세포 점 크기를 조절해서 많은 세포를 보기 쉽게 표현했습니다.

정리하면, 이 코드는 **QC 필터링, 정규화, 고변동 유전자 선택, 스케일링, PCA, Harmony 배치 보정, UMAP, 이웃 탐색, 군집화, 시각화** 순서로 진행됩니다. 그리고 이 코드의 핵심 의의는 여러 샘플의 단일세포 데이터를 단순히 합치는 것이 아니라, **Harmony를 이용해서 샘플 간 기술적 차이를 줄인 뒤 세포 군집 구조를 더 정확하게 해석했다는 점**에 있습니다.
