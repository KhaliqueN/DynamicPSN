from __future__ import division
import collections
from sklearn import metrics
from sklearn.metrics import matthews_corrcoef, precision_recall_fscore_support as score
from sklearn.decomposition import PCA
from collections import Counter
from random import shuffle
from sklearn.linear_model import LogisticRegression
import random
import csv
import time
import math
import pandas as pd
import os
import numpy as np
import warnings
warnings.filterwarnings("ignore")


def read_data(data_directory, file, out_par_idx, NUM_FOLDS_inner, outputDirectory):
    df = pd.read_csv(data_directory + "/" + file, header = None, sep = "\t")

    # Get out_test #
    out_test = pd.read_csv(outputDirectory+"/partitions/prot_idx" + str(out_par_idx+1) + ".txt", header = None, sep = "\t")
    out_test['sort_index'] = list(range(out_test.shape[0]))

    df_out_test = pd.merge(df, out_test, on = [0], how = 'inner')
    df_out_test = df_out_test.sort_values(by = 'sort_index')
    df_out_test = df_out_test.drop(columns = ['sort_index'])

    # Get out_train #
    out_train = pd.read_csv(outputDirectory+"/partitions/tr_prot_idx" + str(out_par_idx+1) + ".txt", header = None, sep = "\t")
    out_train['sort_index'] = list(range(out_train.shape[0]))

    df_out_train = pd.merge(df, out_train, on = [0], how = 'inner')
    df_out_train = df_out_train.sort_values(by = 'sort_index')
    df_out_train = df_out_train.drop(columns = ['sort_index'])

    # Get inner folder # 
    dir1 = outputDirectory+"/partitions/training_partitions/tr_prot_idx" + str(out_par_idx+1) + "/"

    for i in range(NUM_FOLDS_inner):
        # inner test #
        inner_test = pd.read_csv(dir1 + "prot_idx" + str(i+1) + ".txt", header = None, sep = "\t")
        inner_test['sort_index'] = list(range(inner_test.shape[0]))

        globals()['df_inner_test%s' % i] = pd.merge(df_out_train, inner_test, on = [0], how = 'inner')
        globals()['df_inner_test%s' % i] = globals()['df_inner_test%s' % i].sort_values(by = 'sort_index')
        globals()['df_inner_test%s' % i] = globals()['df_inner_test%s' % i].drop(columns = ['sort_index'])

        # inner train #
        inner_train = pd.read_csv(dir1 + "tr_prot_idx" + str(i+1) + ".txt", header = None, sep = "\t")
        inner_train['sort_index'] = list(range(inner_train.shape[0]))

        globals()['df_inner_train%s' % i] = pd.merge(df_out_train, inner_train, on = [0], how = 'inner')
        globals()['df_inner_train%s' % i] = globals()['df_inner_train%s' % i].sort_values(by = 'sort_index')
        globals()['df_inner_train%s' % i] = globals()['df_inner_train%s' % i].drop(columns = ['sort_index'])
    
    return df_out_train, df_out_test, df_inner_test0, df_inner_train0, df_inner_test1, df_inner_train1, df_inner_test2, df_inner_train2, df_inner_test3, df_inner_train3, df_inner_test4, df_inner_train4 




def classifying_LR_l2(data_directory,file,Directory_Save,NUM_FOLDS,outputDirectory):


    NUM_FOLDS = NUM_FOLDS
    # test_indices_outer, train_indices_outer, keys_outer = partition(X, Y, NUM_FOLDS)
    opt_c = [] 


    for k in range(NUM_FOLDS):


        df_out_train, df_out_test, df_inner_test0, df_inner_train0, df_inner_test1, df_inner_train1, df_inner_test2, df_inner_train2, df_inner_test3, df_inner_train3, df_inner_test4, df_inner_train4 = read_data(data_directory, file, k, NUM_FOLDS,outputDirectory)


        # ## testing ######

        # name = Directory_Save + '/tr_' + str(k)+ '.txt'
        # with open(name, 'w') as f:
        #     f.writelines("%s\n" % y for y in np.array(df_out_train.iloc[:, 0]))

        # name = Directory_Save + '/te_' + str(k)+ '.txt'
        # with open(name, 'w') as f:
        #     f.writelines("%s\n" % y for y in np.array(df_out_test.iloc[:, 0]))

        # ###############
        test_outer = np.array(df_out_test.iloc[:, 2:])
        test_outer_labs = np.array(df_out_test.iloc[:, 1])
        gene_test = np.array(df_out_test.iloc[:, 0])

        train_outer = np.array(df_out_train.iloc[:, 2:])
        train_outer_labs = np.array(df_out_train.iloc[:, 1])


        c = np.logspace(-8, 8, num=10, base=2.0)

        c_optimal = 1

        score_max = 0

        score_list_outer = []

        for item in c:
            score_list_inner = []
            for n in range(NUM_FOLDS):

                df_inner_train = globals()['df_inner_train%s' % n]
                df_inner_test = globals()['df_inner_test%s' % n]

                dat_train = np.array(df_inner_train.iloc[:, 2:])
                labs_train = np.array(df_inner_train.iloc[:, 1])
                dat_test = np.array(df_inner_test.iloc[:, 2:])
                labs_test = np.array(df_inner_test.iloc[:, 1])

                model = LogisticRegression(C=item)

                model.fit(dat_train, labs_train)

                # optimizing average accuracy over 5 folds
                score_each_fold = model.score(dat_test, labs_test)
               
                score_list_inner = np.append(score_list_inner, score_each_fold)

            if score_list_inner.mean() > score_max:
                score_max = score_list_inner.mean()
                c_optimal = item


        # training model using optimal hyperparameter ###############
        opt_c = np.append(opt_c, c_optimal)

        model = LogisticRegression(C=c_optimal)
        model.fit(train_outer, train_outer_labs)

        # save accuracy over each fold
        score_each_fold_outer = model.score(test_outer, test_outer_labs)
        score_list_outer = np.append(score_list_outer, score_each_fold_outer)

        y_predict = model.predict(test_outer)
        d = {'Protein':gene_test, 'Actual':test_outer_labs,'Predicted':y_predict}
        predictions = pd.DataFrame(d)

        #save indexes of genes
        name = Directory_Save + '/Prediction' + str(k)+ '.txt'
        predictions.to_csv(name, sep='\t', index=False)



    #save optimal c for each fold
    name = 'Optimal_C.txt'
    with open(Directory_Save + '/' + name, 'w') as f:
        f.writelines("%s\n" % y for y in opt_c)

	# accuracy over each fold 
    name = 'Avg-accuracy.txt'
    with open(Directory_Save + '/' + name, 'w') as f:
        f.writelines("%s\n" % y for y in score_list_outer)



def classify(dataset, NUM_FOLDS, outputDirectory):
    """ main function for doing classification

    Parameters
    ----------
    dataset : str
        Name of dataset
    sampling_type : int
        Type of sampling method
    metric: str
        Metric for evaluation (Accuracy or Matthews correlation coefficient)
    save_acc_perclass: bool
        a flag for whether or not this function saves accuracy per class

    """

    data_directory = dataset
    dict_results = {}
    files = [f for f in os.listdir(data_directory) if f.endswith(".txt")]

    Directory_Save = outputDirectory+"/Classification-results/"
    if not os.path.exists(Directory_Save):
    	os.makedirs(Directory_Save)

    start = time.time()
    classifying_LR_l2(data_directory, files[0], Directory_Save, NUM_FOLDS, outputDirectory)
    end = time.time()
    elapsed = ((end - start) / 60)  # in minutes

    dict_results[files[0]] = [elapsed]

    #save runtime
    # Directory_Save2 = os.path.join(ROOT_PATH, "Results/"+dataset)
    name = 'runtimes_'+os.path.basename(files[0])
    with open(Directory_Save + '/' + name, 'w') as f:
        for key in dict_results:
            f.write(str(key).split('.')[0] + "\t" + str(dict_results[key][0]) + "\n")
