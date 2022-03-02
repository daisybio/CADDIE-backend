from tasks.task_hook import TaskHook


def quick_task(task_hook: TaskHook):
    def run_closeness(parameters):
        from .harmonic_centrality import harmonic_centrality

        def closeness_progress(progress, status):
            task_hook.set_progress(2 / 3 + 1 / 3 * progress, status)

        def closeness_set_result(result):
            task_hook.set_results(result)

        # Prepare intermediate hook
        closeness_task_hook = TaskHook(parameters,
                                       task_hook.data_directory,
                                       closeness_progress,
                                       closeness_set_result)

        # Run closeness centrality
        harmonic_centrality(closeness_task_hook)

    def run_multi_steiner(parameters):
        from .multi_steiner import multi_steiner

        def ms_progress(progress, status):
            task_hook.set_progress(0 + 2 / 3 * progress, status)

        def ms_set_result(result):
            """
            We run closeness on the results of ms tree
            :param result:
            :return:
            """
            # node_attributes = result.get("node_attributes", {})
            cancer_types = result.get("cancer_types", [])
            # node_types = node_attributes.get("node_types", {})
            # seeds are all seed nodes with type see or nodes that were not in node_Types but are still seed nodes
            # TODO how is that possible? do we add nodes but not node_types?
            # seeds = [seed for seed in result["network"]["nodes"] if node_types.get(seed) == 'CancerNode']
            # seeds = result["network"]["nodes"]
            # print('length seeds')
            # print(len(seeds))
            if len(result["network"]["nodes"]) == 0:
                task_hook.set_results({"network": {"nodes": [], "edges": []}})
                return

            print('Multi steiner results')

            closeness_parameters = {
                "seeds": result["network"]["nodes"],
                "result_size": 10,
                "hub_penalty": 1,
                "target": "drugs",
                "include_non_approved_drugs": True,
                "include_indirect_drugs": False,
                'cancer_types': cancer_types,
                'cancer_dataset': result['cancer_dataset'],
                'gene_interaction_dataset': result['gene_interaction_dataset'],
                'drug_interaction_dataset': result['drug_interaction_dataset']
            }
            print("running closeness")
            run_closeness(closeness_parameters)

        parameters["num_trees"] = 1
        parameters["hub_penalty"] = 1

        # Prepare intermediate hook
        ms_task_hook = TaskHook(parameters,
                                task_hook.data_directory,
                                ms_progress,
                                ms_set_result)

        # Run multi_steiner
        print("run multi steiner")
        multi_steiner(ms_task_hook)

    run_multi_steiner(task_hook.parameters)
