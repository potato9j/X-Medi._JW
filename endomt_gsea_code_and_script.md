# EndoMT GSEA 분석 발표 자료

## PPT에 게시할 코드

```r
# 1. 분석 환경 설정 및 데이터 로드
# (필수 라이브러리 로드 및 사전에 처리된 RDS 데이터 파일을 불러옵니다.)
library(Seurat)
library(dplyr)
library(msigdb)
library(fgsea)
library(ggplot2)
library(patchwork)

endo_cells <- readRDS("Glioblastoma_Endo_Before_GSEA.rds")
Idents(endo_cells) <- endo_cells$seurat_clusters

# ---------------------------------------------------------------------
# 2. 발현량 데이터 추출 및 분석 대상 그룹 정의
# ---------------------------------------------------------------------
# 유전자 발현량 매트릭스를 추출하고, 분석 대상(EndoMT)과 대조군(Normal) 세포를 지정합니다.
expr_matrix <- LayerData(endo_cells, assay = "RNA", layer = "data")

# 비교 군집 설정: 실험군(Cluster 12, 13) 및 대조군(Cluster 9)
zombie_cells <- WhichCells(endo_cells, idents = c("12", "13"))
normal_cells <- WhichCells(endo_cells, idents = "9")

# ---------------------------------------------------------------------
# 3. 차등 발현 분석(DEA) 및 유전자 랭킹 산출
# ---------------------------------------------------------------------
# 각 유전자별 Log2 Fold Change(Log2FC)를 직접 산출하여 분석용 랭킹을 생성합니다.
zombie_mean <- log2(rowMeans(expm1(expr_matrix[, zombie_cells])) + 1)
normal_mean <- log2(rowMeans(expm1(expr_matrix[, normal_cells])) + 1)
log2fc <- zombie_mean - normal_mean

# GSEA 분석을 위해 Log2FC 값을 기준으로 유전자를 내림차순 정렬합니다.
ranked_genes <- sort(log2fc, decreasing = TRUE)

# ---------------------------------------------------------------------
# 4. 참조 유전자 세트(Gene Set) 확보 및 GSEA 실행
# ---------------------------------------------------------------------
# MSigDB에서 검증된 Hallmark 및 GO 생물학적 과정(BP) 유전자 세트를 로드합니다.
m_df_h <- msigdb(species = "Homo sapiens", category = "H")
m_df_go <- msigdb(species = "Homo sapiens", category = "C5", subcategory = "BP")

pathways <- list(
  "EMT_Hallmark" = m_df_h %>% filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% pull(gene_symbol),
  "TGF_Beta_Signaling" = m_df_go %>% filter(gs_name == "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY") %>% pull(gene_symbol)
)

# fgsea 알고리즘을 이용하여 유전자 집단 농축 분석을 수행합니다.
set.seed(42)
fgsea_res <- fgsea(pathways = pathways, stats = ranked_genes, minSize = 15, maxSize = 500)

# ---------------------------------------------------------------------
# 5. 분석 결과 시각화 (Enrichment Plot)
# ---------------------------------------------------------------------
# 주요 신호전달 경로의 농축 양상을 시각화하여 형질 전환을 증명합니다.
p1 <- plotEnrichment(pathways[["EMT_Hallmark"]], ranked_genes) +
  labs(title = "1. EMT Hallmark Enrichment", subtitle = "EndoMT 과정의 전형적인 유전자 농축 확인") +
  theme_minimal() + theme(plot.title = element_text(face = "bold", size = 14))

p2 <- plotEnrichment(pathways[["TGF_Beta_Signaling"]], ranked_genes) +
  labs(title = "2. TGF-beta Signaling Enrichment", subtitle = "TGF-beta 신호전달 경로의 활성 양상 분석") +
  theme_minimal() + theme(plot.title = element_text(face = "bold", size = 14))

# 최종 결과 비교 출력
p1 / p2
```

## 발표 대본

이 코드는 혈관 내피세포 하위 군집 데이터인 `endo_cells`를 대상으로, EndoMT와 관련된 분자적 특성이 실제로 유전자 수준에서 농축되어 있는지를 검증하기 위해 GSEA, 즉 유전자 집합 농축 분석을 수행하는 코드입니다. 목적은 EndoMT 성격을 가진 것으로 해석되는 군집과 정상 내피세포 성격의 군집을 비교하여, EMT 관련 유전자 집합과 TGF-beta 신호전달 경로가 실험군에서 더 강하게 나타나는지를 확인하는 것입니다. 먼저 `library(Seurat)`, `library(dplyr)`, `library(msigdb)`, `library(fgsea)`, `library(ggplot2)`, `library(patchwork)`를 불러옵니다. `Seurat`는 단일세포 객체를 다루는 데 사용되고, `dplyr`는 유전자 세트 필터링에 사용되며, `msigdb`는 MSigDB 유전자 세트를 불러오는 역할을 합니다. `fgsea`는 GSEA 알고리즘을 빠르게 수행하기 위한 패키지이고, `ggplot2`와 `patchwork`는 결과 시각화에 사용됩니다. 그다음 `readRDS("Glioblastoma_Endo_Before_GSEA.rds")`를 통해 분석용 RDS 파일을 불러오고, `Idents(endo_cells) <- endo_cells$seurat_clusters`를 이용하여 Seurat 객체의 군집 정보를 현재 분석 기준으로 설정합니다.

이후 `LayerData(endo_cells, assay = "RNA", layer = "data")`를 사용하여 RNA assay의 정규화된 발현량 매트릭스를 추출합니다. 그 다음 `WhichCells()`를 이용해 비교할 두 그룹을 정의합니다. 여기서 `zombie_cells <- WhichCells(endo_cells, idents = c("12", "13"))`는 EndoMT 성격을 가진 실험군으로 설정된 클러스터 12번과 13번 세포를 선택하는 부분이고, `normal_cells <- WhichCells(endo_cells, idents = "9")`는 정상 내피세포 성격의 대조군으로 설정된 클러스터 9번 세포를 선택하는 부분입니다. 즉, 이 코드는 두 집단 간의 발현 차이를 바탕으로 EndoMT와 관련된 경로가 정말로 활성화되어 있는지를 검증하는 구조입니다.

그다음 단계에서는 GSEA에 사용할 유전자 랭킹을 계산합니다. `zombie_mean <- log2(rowMeans(expm1(expr_matrix[, zombie_cells])) + 1)`과 `normal_mean <- log2(rowMeans(expm1(expr_matrix[, normal_cells])) + 1)`는 각각 실험군과 대조군에서 각 유전자의 평균 발현량을 계산하는 부분입니다. 여기서 `expm1()`은 로그 변환된 발현값을 원래 스케일로 되돌린 뒤 평균을 구하고, 다시 `log2(... + 1)`를 적용하여 안정적인 비교가 가능하도록 만든 것입니다. 이후 `log2fc <- zombie_mean - normal_mean`을 통해 각 유전자의 Log2 Fold Change를 직접 계산합니다. 이 값이 양수이면 실험군에서 더 높게 발현된다는 뜻이고, 음수이면 대조군에서 더 높게 발현된다는 뜻입니다. 마지막으로 `ranked_genes <- sort(log2fc, decreasing = TRUE)`를 사용하여 Log2FC 기준으로 유전자를 내림차순 정렬합니다. 이 정렬된 벡터가 바로 GSEA의 입력이 되는 유전자 랭킹입니다.

이후 참조 유전자 세트를 확보하고 실제 GSEA를 수행합니다. `msigdb(species = "Homo sapiens", category = "H")`는 Hallmark 유전자 세트를 불러오고, `msigdb(species = "Homo sapiens", category = "C5", subcategory = "BP")`는 GO Biological Process 유전자 세트를 불러옵니다. 그다음 `pathways` 리스트를 만들어, `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`과 `GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY`에 해당하는 유전자만 각각 추출합니다. 전자는 EMT 관련 표준 유전자 집합이고, 후자는 TGF-beta 수용체 신호전달과 관련된 생물학적 과정 유전자 집합입니다. `set.seed(42)`는 분석 결과의 재현성을 확보하기 위한 설정이고, `fgsea(pathways = pathways, stats = ranked_genes, minSize = 15, maxSize = 500)`는 준비한 유전자 랭킹과 경로 정보를 바탕으로 실제 농축 분석을 수행하는 부분입니다. 여기서 `minSize = 15`와 `maxSize = 500`은 너무 작거나 너무 큰 유전자 세트를 제외하여 해석 가능성을 높이기 위한 조건입니다.

마지막 단계에서는 GSEA 결과를 enrichment plot으로 시각화합니다. `p1 <- plotEnrichment(pathways[["EMT_Hallmark"]], ranked_genes)`는 EMT Hallmark 유전자 세트가 실험군 쪽 상위 랭킹 유전자에 얼마나 집중되어 있는지를 보여주는 그래프이고, `p2 <- plotEnrichment(pathways[["TGF_Beta_Signaling"]], ranked_genes)`는 TGF-beta signaling 관련 유전자 세트의 농축 양상을 보여주는 그래프입니다. 각각 `labs()`를 통해 제목과 부제목을 붙였고, `theme_minimal()`과 `theme(plot.title = element_text(face = "bold", size = 14))`를 사용하여 발표용으로 보기 좋은 형태로 정리했습니다. 마지막으로 `p1 / p2`를 사용하여 두 그래프를 세로로 배치함으로써, EMT 경로와 TGF-beta 경로가 함께 농축되는지를 한 화면에서 비교할 수 있도록 구성했습니다. 정리하면, 이 코드는 EndoMT 성격의 군집과 정상 내피세포 군집 사이의 유전자 발현 차이를 기반으로 GSEA를 수행하여, EMT와 TGF-beta signaling 경로가 실제로 실험군에서 농축되는지를 검증하는 분석 코드입니다. 핵심 의의는 단순히 개별 마커 유전자 발현만 보는 것이 아니라, 관련 유전자 집합 전체의 방향성과 농축 패턴을 통해 EndoMT 해석을 보다 체계적으로 뒷받침한다는 점에 있습니다.
