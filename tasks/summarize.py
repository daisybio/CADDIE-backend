from tasks.task_hook import TaskHook
from caddie.tasks import task_result
from caddie import models
from deepmerge import always_merger
import json
from collections import Counter

def _merge_node_attribute(attribute, task_result_cur, task_result_object):
    for k, v in task_result_cur['node_attributes'][attribute].items():
        if k not in task_result_object['node_attributes'][attribute] or task_result_object['node_attributes'][attribute][k] is None:
            task_result_object['node_attributes'][attribute][k] = v
        else:
            if isinstance(v, bool):
                task_result_object['node_attributes'][attribute][k] = task_result_object['node_attributes'][attribute][k] or v
            elif isinstance(v, str):
                # there are just 'Drug', 'Node' and 'CancerNode'. Drugs are always Drug, if a node has been CancerNode once, set it to 'CancerNode'
                task_result_object['node_attributes'][attribute][k] = v if v == task_result_object['node_attributes'][attribute][k] else 'CancerNode'
            elif isinstance(v, float):
                task_result_object['node_attributes'][attribute][k] = (task_result_object['node_attributes'][attribute][k] + v) / 2
            elif isinstance(v, int):
                task_result_object['node_attributes'][attribute][k] = (task_result_object['node_attributes'][attribute][k] + v) / 2
    return task_result_object


def summarize(task_hook: TaskHook):
    """create summary network and summary statistics

    Args:
        task_hook (TaskHook): Task object for summarization
    """
    source_tasks = task_hook.parameters["source_tasks"]
    task_result_object = task_result(models.Task.objects.get(token=source_tasks.pop()))
    for token_str in source_tasks:
        task = models.Task.objects.get(token=token_str)
        task_result_cur = task_result(task)
        # values are needed to be merged separately, always_merger only takes first argument
        task_result_object = _merge_node_attribute('node_types', task_result_cur, task_result_object)
        task_result_object = _merge_node_attribute('degrees', task_result_cur, task_result_object)
        task_result_object = _merge_node_attribute('is_seed', task_result_cur, task_result_object)
        task_result_object = _merge_node_attribute('is_result', task_result_cur, task_result_object)

        if 'scores' in task_result_cur:
            task_result_object = _merge_node_attribute('scores', task_result_cur, task_result_object)

        node_attributes = task_result_object['node_attributes']
        del task_result_cur['node_attributes']
        del task_result_object['node_attributes']
        task_result_object = always_merger.merge(task_result_object, task_result_cur)
        task_result_object['node_attributes'] = node_attributes

    # fetch node names to avoid repetetive database lookups
    node_names = {}
    for node in task_result_object['network']['nodes']:
        if node[0] == 'd':
            node_names[node] = models.Drug.objects.get(id=node[1:]).name
        else:
            node_names[node] = models.Gene.objects.get(id=node[1:]).name

    # count drug occurences
    task_result_object['drug_counts'] = dict(Counter([node_names[node] for node in task_result_object['network']['nodes']if node[0]=='d']))
    # count gene occurences of non-seed genes
    node_counts_counter = Counter([node for node in task_result_object['network']['nodes']]) #  if not task_result_object['node_attributes']['is_seed'][node]
    gene_counts = {node_names[node]: count for node, count in node_counts_counter.items() if task_result_object['node_attributes']['node_types'][node] == 'Node'}
    cancer_gene_counts = {node_names[node]: count for node, count in node_counts_counter.items() if task_result_object['node_attributes']['node_types'][node] == 'CancerNode'}
    task_result_object['gene_counts'] = gene_counts
    task_result_object['cancer_gene_counts'] = cancer_gene_counts
    
    traces_degree = {'Drug': {'x': [], 'y': [], 'names': []}, 'CancerNode': {'x': [], 'y': [], 'names': []}, 'Node': {'x': [], 'y': [], 'names': []}}
    if 'scores' in task_result_object['node_attributes']:
        for node in task_result_object['network']['nodes']:
            degree = task_result_object['node_attributes']['degrees'][node]
            score = task_result_object['node_attributes']['scores'][node] if node in task_result_object['node_attributes']['scores'] else None
            node_type = task_result_object['node_attributes']['node_types'][node]
            traces_degree[node_type]['x'].append(degree)
            traces_degree[node_type]['y'].append(score)
            traces_degree[node_type]['names'].append(node_names[node])
    task_result_object['traces_degree'] = traces_degree

    # remove duplicate nodes
    task_result_object['network']['nodes'] = list(set(task_result_object['network']['nodes']))
    task_hook.set_results(task_result_object)
