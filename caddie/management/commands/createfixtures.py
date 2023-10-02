import os

import django
from django.core.management import BaseCommand
import subprocess

django.setup()


class Command(BaseCommand):
    help = 'Creates fixtures'

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        
        su_username = os.environ.get('DJANGO_SUPERUSER_USERNAME', 'nope')
        su_email = os.environ.get('DJANGO_SUPERUSER_EMAIL', '')
        su_password = os.environ.get('DJANGO_SUPERUSER_PASS', 'nope123')
        
        process = subprocess.Popen(f'python3 manage.py createsuperuser --noinput --email {su_email} --username {su_username}', shell=True, stdout=subprocess.PIPE)
        process.wait()
        print('RETURNCODE')
        print(process.returncode)


        print(su_username, su_password)
        
        # from django.contrib.auth.models import User
        # user, created = User.objects.get_or_create(username=su_username)
        # user.set_password(su_password)
        # user.is_superuser = True
        # user.is_staff = True
        # user.save()

        # if created:
        #     print('created caddie superuser')
