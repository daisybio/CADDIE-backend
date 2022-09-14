from django.contrib import admin
import caddie.models as models

admin.site.register(models.Gene)
admin.site.register(models.Drug)
admin.site.register(models.CancerDataset)
admin.site.register(models.GeneGeneInteraction)
admin.site.register(models.CancerType)
admin.site.register(models.CancerGeneEntity)
admin.site.register(models.GeneDrugInteraction)
admin.site.register(models.Task)
admin.site.register(models.DrugDataset)
admin.site.register(models.GeneExpressionLevel)
admin.site.register(models.MutationCancerType)
