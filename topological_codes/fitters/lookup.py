def lookuptable_decoding(training_results, real_results):
    """
    Calculates the logical error probability using postselection decoding.
    This postselects all results with trivial syndrome.
    Args:
        training_results (dict): A results dictionary, as produced by the
            ``process_results`` method of a code.
        real_results (dict): A results dictionary, as produced by the
            ``process_results`` method of a code.
    Returns:
        dict: Dictionary of logical error probabilities for
            each of the encoded logical states whose results were given in
            the input.
    Additional information:
        Given a two dictionaries of results, as produced by a code object,
        thelogical error probability is calculated for lookup table
        decoding. This is done using `training_results` as a guide to which
        syndrome is most probable for each logical value, and the
        probability is calculated for the results in `real_results`.
    """

    logical_prob = {}
    for log in real_results:
        shots = 0
        incorrect_shots = 0
        for string in real_results[log]:

            p = {}
            for testlog in ["0", "1"]:
                if string in training_results[testlog]:
                    p[testlog] = training_results[testlog][string]
                else:
                    p[testlog] = 0

            shots += real_results[log][string]
            if p["1" * (log == "0") + "0" * (log == "1")] > p[log]:
                incorrect_shots += real_results[log][string]

        logical_prob[log] = incorrect_shots / shots

    return logical_prob
