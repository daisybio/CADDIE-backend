import random
import time

from tasks.task_hook import TaskHook


def sample_task(task_hook: TaskHook):
    file_path = task_hook.data_directory + 'test_data.csv'
    seeds = task_hook.seeds

    for i in range(20):
        task_hook.set_progress(i / 20.0, 'In the loop')
        time.sleep(random.randint(1, 10) / 100.0)

    task_hook.set_results({
        'network': {'nodes': ['Q9H4P4', 'P00533'], 'edges': [{'from': 'Q9H4P4', 'to': 'P00533'}]},
    })
