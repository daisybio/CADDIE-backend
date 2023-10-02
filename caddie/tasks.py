import json
from datetime import datetime
import caddie.settings as settings
import redis
import rq
import os

from tasks.task_hook import TaskHook

import django
django.setup()

qr_r = redis.Redis(**settings.redis_config)
rq_tasks = rq.Queue('caddie_tasks', connection=qr_r)

r = redis.Redis(**settings.redis_config_decode)


def run_task(token, algorithm, parameters_json):
    def set_progress(progress, status):
        r.set(f'{token}_progress', f'{progress}')
        r.set(f'{token}_status', f'{status}')

    def set_result(results):
        r.set(f'{token}_result', json.dumps(results, allow_nan=True))
        r.set(f'{token}_finished_at', f'{datetime.now().timestamp()}')
        r.set(f'{token}_done', '1')

        set_progress(1.0, 'Done.')

    set_progress(0.0, 'Computation started')

    worker_id = os.getenv('RQ_WORKER_ID')
    r.set(f'{token}_worker_id', f'{worker_id}')
    job_id = os.getenv('RQ_JOB_ID')
    r.set(f'{token}_job_id', f'{job_id}')
    r.set(f'{token}_started_at', f'{datetime.now().timestamp()}')

    parameters = json.loads(parameters_json)
    task_hook = TaskHook(parameters, 'data/networks/', set_progress, set_result)
    try:
        if algorithm == 'dummy':
            raise RuntimeError('Dummy algorithm for testing purposes.')
        elif algorithm in ['multisteiner', 'exampledrugtarget', 'multi-steiner']:
            from tasks.multi_steiner import multi_steiner
            multi_steiner(task_hook)
        elif algorithm == 'keypathwayminer':
            from tasks.keypathwayminer_task import kpm_task
            kpm_task(task_hook)
        elif algorithm in ['trustrank', 'exampledrug', 'simpledrug']:
            from tasks.trust_rank import trust_rank
            trust_rank(task_hook)
        elif algorithm in ['harmonic', 'harmonic_centrality']:
            from tasks.harmonic_centrality import harmonic_centrality
            harmonic_centrality(task_hook)
        elif algorithm in ['degree', 'degree_centrality']:
            from tasks.degree_centrality import degree_centrality
            degree_centrality(task_hook)
        elif algorithm in ['proximity', 'network_proximity']:
            from tasks.network_proximity import network_proximity
            network_proximity(task_hook)
        # elif algorithm == 'domino':
        #     from tasks.domino import domino_task
        #     domino_task(task_hook, token)
        elif algorithm in ['betweenness', 'betweenness_centrality']:
            from tasks.betweenness_centrality import betweenness_centrality
            betweenness_centrality(task_hook)
        elif algorithm in ['quick', 'super']:
            from tasks.quick_task import quick_task
            quick_task(task_hook)
        elif algorithm == 'summary':
            from tasks.summarize import summarize
            summarize(task_hook)
    except Exception as e:
        r.set(f'{token}_status', f'{e}')
        r.set(f'{token}_failed', '1')


def refresh_from_redis(task):
    task.worker_id = r.get(f'{task.token}_worker_id')
    if not task.worker_id:
        return

    task.job_id = r.get(f'{task.token}_job_id')
    task.progress = float(r.get(f'{task.token}_progress'))
    task.done = True if r.get(f'{task.token}_done') else False
    task.failed = True if r.get(f'{task.token}_failed') else False
    status = r.get(f'{task.token}_status')
    if not status or len(status) < 255:
        task.status = status
    else:
        task.status = status[:255]
    started_at = r.get(f'{task.token}_started_at')
    if started_at:
        task.started_at = datetime.fromtimestamp(float(started_at))
    finished_at = r.get(f'{task.token}_finished_at')
    if finished_at:
        task.finished_at = datetime.fromtimestamp(float(finished_at))
    task.result = r.get(f'{task.token}_result')


def start_task(task):
    job = rq_tasks.enqueue(run_task, task.token, task.algorithm, task.parameters, job_timeout=30*60*5)
    task.job_id = job.id


def task_stats(task):
    pos = 1
    for j in rq_tasks.jobs:
        if j.id == task.job_id:
            break
        pos += 1

    return {
        'queueLength': rq_tasks.count,
        'queuePosition': pos,
    }


def task_result(task):
    if not task.done:
        return None
    return json.loads(task.result, parse_constant=lambda c: None)
