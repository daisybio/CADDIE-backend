import os

import django
from django.core.management import BaseCommand

django.setup()


class Command(BaseCommand):
    help = 'Creates fixtures'

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        su_username = os.environ.get('SUPERUSER_NAME', 'admin')
        su_password = os.environ.get('SUPERUSER_PASS', 'admin')

        from django.contrib.auth.models import User
        user, created = User.objects.get_or_create(username=su_username)
        user.set_password(su_password)
        user.is_superuser = True
        user.is_staff = True
        user.save()

        if created:
            print('created admin superuser')
