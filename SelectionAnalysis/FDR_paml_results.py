__author__ = 'juan'
   #Perform FDR, based on Benjamini approach, at 5%. Add this to the qvalue in the dictionary
    #Based on this post:
    #http://stats.stackexchange.com/questions/870/multiple-hypothesis-testing-correction-with-benjamini-hochberg-p-values-or-q-va

    total_tests = len(cluster_paml_results)  # Total number of performed tests
    position = 1
    prev_adjusted_pvalue = 0

    print cluster_paml_results
    print total_tests

    sorted_results = sorted(cluster_paml_results, key=itemgetter(4))


    for entry in sorted_results:
        adjusted_pvalue = float(entry[4]) * (total_tests / position)

        print adjusted_pvalue

        #If the value is greater than 1, we set as one (0 < p < 1)
        adjusted_pvalue = min(adjusted_pvalue, 1)

        #Check that the value is not greater than the previous one
        adjusted_pvalue = max(adjusted_pvalue, prev_adjusted_pvalue)

        prev_adjusted_pvalue = adjusted_pvalue
        position += 1

        entry[5] = adjusted_pvalue



    print cluster_paml_results

    #print cluster_paml_results
    #Summary information