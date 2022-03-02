from django.utils import timezone
from django.core.management.base import BaseCommand

from caddie.models import Task


class Command(BaseCommand):
    def add_arguments(self, parser):
        pass

    def handle(self, *args, **kwargs):
        print('Task cleanup...\n')

        print(f'Closing {Task.objects.filter(started_at__isnull=True).count()} queued tasks')
        Task.objects.filter(started_at__isnull=True).update(started_at=timezone.now(), finished_at=timezone.now(),
                                                            failed=True, status='Removed due to server redeployment.')

        print(f'Closing {Task.objects.filter(finished_at__isnull=True).count()} processing tasks')
        Task.objects.filter(finished_at__isnull=True).update(finished_at=timezone.now(), failed=True,
                                                             status='Stopped due to server redeployment.')
