def calculate_metrics(ground_truth_genes: set, amr_genes: set) -> tuple[float, float, float]:
    """
    Retrieves the metrics from the amr and ground truth gene sets.
    :param ground_truth_genes: Set of ground truth genes
    :param amr_genes: Set of amr genes to compare
    :return: four metrics: precision, recall, f1
    """
    true_positives = amr_genes.intersection(ground_truth_genes)
    false_positives = amr_genes.difference(ground_truth_genes)
    false_negatives = ground_truth_genes.difference(amr_genes)
    try:
        precision = len(true_positives) / (len(true_positives) + len(false_positives))
    except ZeroDivisionError:
        precision = 0.0
    try:
        recall = len(true_positives) / (len(true_positives) + len(false_negatives))
    except ZeroDivisionError:
        recall = 0.0
    try:
        f1_score = 2 / (precision ** -1 + recall ** -1)
    except ZeroDivisionError:
        f1_score = 0.0
    return precision, recall, f1_score
