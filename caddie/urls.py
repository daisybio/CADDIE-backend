"""caddie URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/3.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include'blog.urls'))
"""


from django.contrib import admin
from django.urls import path, include

import caddie.views as views
from django.views.decorators.cache import cache_page

# cache time is 6 hours
urlpatterns = [
    path('network/', views.Network.as_view()),
    path('cancer_datasets/', cache_page(21600)(views.CancerDatasetView.as_view())),
    path('interaction_gene_datasets/', cache_page(21600)(views.InteractionGeneDatasetView.as_view())),
    path('interaction_drug_datasets/', cache_page(21600)(views.InteractionDrugDatasetView.as_view())),
    path('cancer_types/', cache_page(21600)(views.CancerTypeView.as_view())),
    path('related_cancer_types/', views.RelatedCancerTypeView.as_view()),
    path('comorbidities_node/', views.ComorbidityNodeView.as_view()),
    path('comorbidities_cancer_type/', views.ComorbidityCancerTypeView.as_view()),
    path('task_summarize/', views.task_summarize),
    path('task/', views.TaskView.as_view()),
    path('tasks/', views.tasks_view),
    path('task_result/', views.result_view),
    path('summary/', views.SummaryView.as_view()),
    path('nodeRelations/', views.NodeInteractionLookup.as_view()),
    path('graph_export/', views.graph_export),
    path('diseases/', views.DiseaseView.as_view()),
    path('query_disease_genes/', views.DiseaseGeneInteractionView.as_view()),
    path('query_nodes/', views.query_nodes),
    path('query_tissue_genes/', views.query_tissue_genes),
    path('query_expression_cancer_type_genes/', views.query_expression_cancer_type_genes),
    path('query_mutation_cancer_type_genes/', views.query_mutation_cancer_type_genes),
    path('drug_interactions/', cache_page(21600)(views.GeneDrugInteractionView.as_view())),
    path('mutation_scores/', cache_page(21600)(views.MutationScoreView.as_view())),
    path('tissue_expression/', cache_page(21600)(views.TissueExpressionView.as_view())),
    path('tissues/', cache_page(21600)(views.TissueView.as_view())),
    path('gene_expression/', cache_page(21600)(views.GeneExpressionView.as_view())),
    path('drug_target_actions/', (views.DrugTargetActionView.as_view())),
    path('expression_cancer_types/', cache_page(21600)(views.ExpressionCancerTypeView.as_view())),
    path('mutation_cancer_types/', cache_page(21600)(views.MutationCancerTypeView.as_view())),
    path('drug_status/', cache_page(21600)(views.DrugStatusView.as_view())),
    path('drug_interaction_lookup/', views.DrugInteractionLookup.as_view()),
    path('cancernet_lookup/', views.CancernetLookupView.as_view()),
    path('gene_drug_interaction_lookup/', views.GeneDrugInteractionLookup.as_view()),
    path('vcf_lookup/', views.VCFSeedLookup.as_view()),
    path('admin/clearcache/', include('clearcache.urls')),
    path('admin/', admin.site.urls),

    path('cancer_genes_internal/', views.InternalCancerGeneView.as_view()),
    path('cancer_gene_interactions_internal/', views.InternalGeneInteractionView.as_view()),
    path('cancer_drug_interactions_internal/', views.InternalGeneDrugInteractionView.as_view()),
]
