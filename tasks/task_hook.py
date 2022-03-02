class TaskHook:

    def __init__(self, parameters, data_directory, set_progress, set_result):
        self.__parameters = parameters
        self.__data_directory = data_directory
        self.__set_progress = set_progress
        self.__set_result = set_result

    @property
    def seeds(self):
        """
        Returns seeds selected for the algorithm.

        :return: Proteins as a list of Uniprot accession codes (e.g. ["P61970", "Q9H4P4"])
        """
        return self.__parameters['seeds']

    @property
    def parameters(self):
        """
        Returns parameters selected for the algorithm.

        :return: Parameters as dictionary (e.g. {"proteins": [...], "paramA": 123, "paramB": True, ...})
        """
        return self.__parameters

    @property
    def data_directory(self):
        """
        Returns the data directory including trailing slash.

        :return: Data directory (e.g. '/app/data/')
        """
        return self.__data_directory

    def set_progress(self, progress, status):
        """
        To be called to indicate computation progress.

        :param progress: A float between 0.0 and 1.0 (e.g. 0.5 for 50% progress)
        :param status: A string indicating the status (e.g. 'Parsing file')
        :return:
        """
        self.__set_progress(progress, status)

    def set_results(self, results):
        """
        To be called when the computation is finished.

        :param results: A dictionary containing a networks entry, each network having nodes and edges.
        (e.g. {"network": {"nodes": ["P61970", "Q9H4P4"], "edges": [{"from": "P61970", "to": "Q9H4P4"}]}})
        """
        self.__set_result(results)
