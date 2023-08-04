from keras.models import Sequential, Model, load_model
from keras.layers import Dense, Dropout, Activation, Input, Lambda, LeakyReLU, Concatenate, BatchNormalization
from keras.optimizers import Adam
import keras.backend as K
import pandas as pd
import os
import numpy as np
from sklearn import metrics
import random
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import argparse as ap
from keras.utils import np_utils
import math
import config

class clf():
    
    def __init__(self, n_in=6523, hidden=[1200, 500, 300, 200], drop_out=[0, 0, 0, 0, 0], batch_size=512, kernel_initializer='glorot_normal', learning_rate=0.001, splicing_amount=True, CIRIdeepA=False):
                
        if splicing_amount:
            n_in += 2
        input = Input(shape=(n_in, ), dtype='float32', name='input') 
        firstlayer = True
        if drop_out[0] > 0: 
            h = Dropout(drop_out[0])(input)
            firstlayer = False
        for i in range(len(hidden)):
            if firstlayer:
                h = Dense(hidden[i], activation=None, kernel_initializer=kernel_initializer)(input) 
                h = BatchNormalization()(h)
                h = Activation('relu')(h)
                firstlayer = False
            else:
                h = Dense(hidden[i], activation=None, kernel_initializer=kernel_initializer)(h)
                h = BatchNormalization()(h)
                h = Activation('relu')(h)
            if drop_out[i+1] > 0:
                h = Dropout(drop_out[i+1])(h)
        if not CIRIdeepA:
            output = Dense(1, activation='sigmoid', name='output')(h)
            model = Model(inputs=input, outputs=output)
            adam = Adam(lr=learning_rate)
            model.compile(optimizer=adam, loss='binary_crossentropy', metrics=['accuracy'])
        else:
            output = Dense(3, activation='softmax', name='output')(h)
            model = Model(inputs=input, outputs=output)
            adam = Adam(lr=learning_rate)
            model.compile(optimizer=adam, loss='categorical_crossentropy', metrics=['accuracy'])

        self.model = model
        self.batch_size = batch_size
        self.n_in = n_in
        self.hidden = hidden
        self.drop_out = drop_out
        self.learning_rate = learning_rate
        self.CIRIdeepA = CIRIdeepA
        
    def fit(self, inputdata):
    
        best_metrics = train_on_balanced_batch(self.model, inputdata, batch_size=self.batch_size, CIRIdeepA=self.CIRIdeepA)
    
        return best_metrics
        
    def predict(self, test_data):
    
        Y_pred = self.model.predict([test_data['X']])
        
        return Y_pred
    
def generate_data_batch(train_list, test_data, seqFeature_df, geneExp_absmax, geneExp_colnames, odir, RBP_dir, splicing_max='', splicing_dir='', splicing_amount=False, n_epoch=100, CIRIdeepA=False):
    '''
    generate dataset for training while fill test data
    '''
    
    # epoch looping
    
    # test_data_saved = False
    n_epoch = config.n_epoch
    
    for i in range(n_epoch):
        
        sample_list = np.array(list(train_list.index))
        
        # 1 epoch
        while len(sample_list) > 0:
            
            # randomly select 4 sample comparisons
            random.seed(config.random_seed_eval)
            idx_inbatch = random.sample(range(len(sample_list)), min(4, len(sample_list)))
            sample_inbatch = sample_list[idx_inbatch]
            sample_list = np.delete(sample_list, idx_inbatch)
            print('epoch: %i; sample in batch: %s' %(i, ', '.join(sample_inbatch)))
            
            X_train_list = []
            Y_train_list = []
            rownames_train_list = []
            X_test_list = []
            Y_test_list = []
            rownames_test_list = []
            
            for sample in sample_inbatch:
            
                label_fn = train_list.loc[sample].label_fn
                X_tmp, Y_tmp, rownames_tmp, colnames = \
                    construct_data_from_label(label_fn, sample, seqFeature_df, geneExp_absmax, geneExp_colnames, RBP_dir, splicing_dir=splicing_dir, splicing_max=splicing_max, splicing_amount=False, phase='train', CIRIdeepA=CIRIdeepA)
                inputdata = {'X':X_tmp, 'Y':Y_tmp, 'rownames':rownames_tmp}
                X_train_tmp, Y_train_tmp, X_test_tmp, Y_test_tmp, rownames_train_tmp, rownames_test_tmp = \
                    split_testdata_out(inputdata, CIRIdeepA=CIRIdeepA)
                    
                X_train_list.append(X_train_tmp)
                Y_train_list.append(Y_train_tmp)
                rownames_train_list.append(rownames_train_tmp)
                X_test_list.append(X_test_tmp)
                Y_test_list.append(Y_test_tmp)
                rownames_test_list.append(rownames_test_tmp)

            # concatenate 4 data
            X_train = np.vstack(X_train_list)
            if not CIRIdeepA:
                Y_train = np.hstack(Y_train_list)
                Y_test = np.hstack(Y_test_list)
            else:
                Y_train = np.vstack(Y_train_list)
                Y_test = np.vstack(Y_test_list)
            rownames_train = np.hstack(rownames_train_list)
            X_test = np.vstack(X_test_list)
            rownames_test = np.hstack(rownames_test_list)

            # add test dataset into test_data in epoch 0
            if i == 0:
                test_data['X'].append(X_test)
                test_data['Y'].append(Y_test)
                test_data['rownames'].append(rownames_test)
            
            # generate train dataset for 4 samples
            train_data = {'X':X_train, 'Y':Y_train, 'rownames':rownames_train, 'colnames':colnames}
            yield train_data

    return
            

def split_testdata_out(inputdata, circid_test=[], test_prop=0.05, random_seed=12345, CIRIdeepA=False):
    '''
    split data into test data and train data
    '''
    test_prop = config.test_prop
    random_seed = config.random_seed_test
    
    X = inputdata['X']
    Y = inputdata['Y']
    rownames = inputdata['rownames']
    
    if not CIRIdeepA:
        # half pos half neg in test data
        test_size = int(len(rownames) * test_prop * 0.5) 
        pos_idx = np.where(Y == 1)[0]
        neg_idx = np.where(Y == 0)[0]
        
        # randomly select test data
        np.random.seed(random_seed)
        np.random.shuffle(pos_idx) ### random seed set
        np.random.seed(random_seed)
        np.random.shuffle(neg_idx) ### random seed set
        
        pos_idx_test = pos_idx[:test_size]
        pos_idx_train = pos_idx[test_size:]
        neg_idx_test = neg_idx[:test_size]
        neg_idx_train = neg_idx[test_size:]
        
        X_test = np.vstack((X[pos_idx_test, ], X[neg_idx_test, ]))
        X_train = np.vstack((X[pos_idx_train, ], X[neg_idx_train, ]))
        Y_test = np.asarray([1.]*len(pos_idx_test) + [0.]*len(neg_idx_test))
        Y_train = np.asarray([1.]*len(pos_idx_train) + [0.]*len(neg_idx_train))
        rownames_test = np.hstack((rownames[pos_idx_test], rownames[neg_idx_test]))
        rownames_train = np.hstack((rownames[pos_idx_train], rownames[neg_idx_train]))
   
    else:
        test_size = int(len(rownames) * test_prop)
        np.random.seed(random_seed)
        idx_shuffled = list(range(len(rownames)))
        for i in range(5):
            np.random.shuffle(idx_shuffled)
        test_idx = idx_shuffled[:test_size]
        train_idx = idx_shuffled[test_size:]
        X_test = X[test_idx, ]
        X_train = X[train_idx, ]
        Y_test = Y[test_idx, ]
        Y_train = Y[train_idx, ]
        rownames_test = rownames[test_idx]
        rownames_train = rownames[train_idx]
    
    return X_train, Y_train, X_test, Y_test, rownames_train, rownames_test


def split_eval_train_data(inputdata, n_fold=5): #######
    '''
    split train data into train and validation
    '''
    X = inputdata['X']
    Y = inputdata['Y']
    rownames = inputdata['rownames']
    
    # randomly select 1/5 data from whole dataset
    fold_size = int(len(rownames)/n_fold)
    idx = np.asarray(range(0, len(rownames)))
    np.random.seed(config.random_seed_eval)
    np.random.shuffle(idx)
    
    idx_val = idx[:fold_size]
    idx_train = idx[fold_size:]
    
    X_val = X[idx_val, ]
    X_train = X[idx_train, ]
    Y_val = Y[idx_val, ]
    Y_train = Y[idx_train, ]
    
    rownames_val = rownames[idx_val]
    rownames_train = rownames[idx_train]
    
    return X_train, X_val, Y_train, Y_val, rownames_train, rownames_val
    

def split_data_into_balanced_minibatch(Y_train, batch_size=64, pos_prop=0.5, CIRIdeepA=False):
    '''
    split train data into balanced mini-batches
    '''
    pos_prop = config.pos_prop
    
    if not CIRIdeepA:
        pos_idx = np.where(Y_train == 1)[0] # shuffle
        neg_idx = np.where(Y_train == 0)[0]
    else:
        neg_idx = np.where(Y_train[:,0] == 1)[0]
        pos_idx = np.where(Y_train[:,0] == 0)[0]
    
    if len(pos_idx) + len(neg_idx) < batch_size:
        batch_size = 2 * min(len(pos_idx), len(neg_idx))
    
    # half pos half neg in batch
    pos_size = int(batch_size * pos_prop)
    neg_size = batch_size - pos_size
    
    pos_batch = len(pos_idx)//pos_size
    neg_batch = len(neg_idx)//neg_size
    n_batch = max(pos_batch, neg_batch)
    
    for i in range(5):
        np.random.shuffle(pos_idx)
    pos_idx = np.tile(pos_idx, n_batch)
    
    for i in range(5):
        np.random.shuffle(neg_idx)
    neg_idx = np.tile(neg_idx, n_batch)
    
    idx_inbatch_list = []
    for i in range(0, n_batch):
        pos_idx_inbatch = pos_idx[i*pos_size:(i+1)*pos_size]
        neg_idx_inbatch = neg_idx[i*neg_size:(i+1)*neg_size]
        idx_inbatch = np.hstack((pos_idx_inbatch, neg_idx_inbatch))
        np.random.shuffle(idx_inbatch)
        idx_inbatch_list.append(idx_inbatch)
        
    return idx_inbatch_list

###################################
def train_on_balanced_batch(model, inputdata, batch_size=64, validation_freq=10, CIRIdeepA=False):
    '''
    train and validation (for each sample pair)
    '''
    
    # construct balanced minibatch
    Y_train = inputdata['Y_train']
    idx_inbatch_list = split_data_into_balanced_minibatch(Y_train, batch_size=batch_size, CIRIdeepA=CIRIdeepA)
    
    # initiate
    best_roc = 0
    best_loss = float('Inf')
    n_batch = len(idx_inbatch_list)
    
    # train on batch for n times
    for i in range(0, n_batch): 
        
        X_train, Y_train, X_val, Y_val, rownames_train, rownames_val = \
            inputdata['X_train'], inputdata['Y_train'], inputdata['X_val'], inputdata['Y_val'], inputdata['rownames_train'], inputdata['rownames_val']
        loss_on_batch = model.train_on_batch(X_train[idx_inbatch_list[i], ], Y_train[idx_inbatch_list[i], ]) # debug: loss unstable
        print('loss on batch%i: %.3f' % (i, loss_on_batch[0]))
        print('accuracy on batch%i: %.3f' % (i, loss_on_batch[1]))
            
        if i == n_batch - 1:
            
            Y_predict = model.predict(X_val)
            eval_val = model.evaluate(X_val, Y_val, verbose=1) 
            eval_train = model.evaluate(X_train, Y_train, verbose=1) 
            
            if not CIRIdeepA:
                roc = metrics.roc_auc_score(Y_val, Y_predict)
                pr = metrics.average_precision_score(Y_val, Y_predict)
            else:
                Y_predict_categorical = (Y_predict == Y_predict.max(axis=1, keepdims=True)).astype('float')
                precision_score_macro = metrics.precision_score(Y_val, Y_predict_categorical, average="macro")
                precision_score_weighted = metrics.precision_score(Y_val, Y_predict_categorical, average="weighted")
                recall_macro = metrics.recall_score(Y_val, Y_predict_categorical, average='macro')
                recall_weighted = metrics.recall_score(Y_val, Y_predict_categorical, average="weighted")
                F1_macro = metrics.f1_score(Y_val, Y_predict_categorical, average='macro')
                F1_weighted = metrics.f1_score(Y_val, Y_predict_categorical, average='weighted')

    print('*' * 50)
    if not CIRIdeepA:
        print('batch %i: loss_val - %.3f; accuracy_val - %.3f; loss_train - %.3f; accuracy_train - %.3f; auroc_val - %.3f; aupr_val - %.3f' % (i, eval_val[0], eval_val[1], eval_train[0], eval_train[1], roc, pr))
    else:
        print('batch %i: loss_val - %.3f; accuracy_val - %.3f; loss_train - %.3f; accuracy_train - %.3f' % (i, eval_val[0], eval_val[1], eval_train[0], eval_train[1]))
    print('*' * 50)
    
    if not CIRIdeepA:
        return {'roc':roc, 'loss_eval':eval_val, 'loss_train':eval_train, 'aupr':pr}
    else:
        return {'loss_train':eval_train[0], 'acc_train':eval_train[1], 'loss_val':eval_val[0], 'acc_val':eval_val[1], 
        'precision_score_macro':precision_score_macro, 'precision_score_weighted':precision_score_weighted, 
        'recall_macro':recall_macro, 'recall_weighted':recall_weighted, 
        'F1_macro':F1_macro, 'F1_weighted':F1_weighted}

###################################
def read_label_fn(label_fn, min_read_cov=20, significance=0.1, CIRIdeepA=False):
    '''
    get circid and label for one sample comparison
    return pos(dict), neg(dict)
    '''
    if not CIRIdeepA:
        
        pos = {}
        neg = {}
        
        with open(label_fn, 'r') as f:
        
            firstline=True
            
            for line in f:
            
                ele = line.strip().split('\t')
                
                if firstline:
                
                    header = {ele[i]:i for i in range(len(ele))} 
                    firstline = False
                    continue
                
                post_pr = float(ele[header['post_pr']])
                eid = ele[header['ID']]
                I1 = int(ele[header['I1']])
                I2 = int(ele[header['I2']])
                S1 = int(ele[header['S1']])
                S2 = int(ele[header['S2']])
                delta = float(ele[header['delta.mle']])
                
                    
                if I1+S1 <= min_read_cov or I2+S2 <= min_read_cov:
                    continue
                if post_pr > 1-significance and min(S1, S2) > 2:
                    pos[eid] = post_pr
                elif post_pr < significance and min(I1, S1)>2 and min(I2, S2)>2:
                    neg[eid] = post_pr
                    
        return pos, neg
    
    else:
        no_diff = {}
        a_higher = {}
        b_higher = {}
        
        with open(label_fn, 'r') as f:
        
            firstline=True
            
            for line in f:
            
                ele = line.strip().split('\t')
                
                if firstline:
                
                    header = {ele[i]:i for i in range(len(ele))} 
                    firstline = False
                    continue
                
                post_pr = float(ele[header['post_pr']])
                eid = ele[header['ID']]
                I1 = int(ele[header['I1']])
                I2 = int(ele[header['I2']])
                S1 = int(ele[header['S1']])
                S2 = int(ele[header['S2']])
                
                #***direction***#
                delta = float(ele[header['delta.mle']])
                
                if I1+S1 <= min_read_cov or I2+S2 <= min_read_cov:
                    continue
                
                if post_pr > 1-significance and min(S1, S2) > 2 and delta < 0:
                    a_higher[eid] = post_pr
                
                if post_pr > 1-significance and min(S1, S2) > 2 and delta > 0:
                    b_higher[eid] = post_pr
                
                elif post_pr < significance and min(I1, S1)>2 and min(I2, S2)>2:
                    no_diff[eid] = post_pr

        return no_diff, a_higher, b_higher

def read_label_fn_predict(label_fn):
    with open(label_fn, 'r') as f:
        circlist = f.read().strip().split('\n')
    return circlist

def read_train_list(fn):
    '''
    return dataframe(sampleID label_fn)
    '''
    df = pd.read_table(fn, index_col=0)    # ID label_fn num_pos num_neg geneExp_fn
    return df

def read_sequence_feature(seqFeature_fn):

    data = pd.read_csv(seqFeature_fn, sep='\t', index_col=0)
    data = data.iloc[:, 13:]

    return data
    
def read_geneExp_absmax(fn):

    df = pd.read_csv(fn, index_col=0, sep='\t')['max']
    colnames = list(df.index)
    rbp_max = df.values

    return rbp_max, colnames

def read_geneExp(sample, geneExp_absmax, nrow=1, phase='train', RBPexp_dir=''):

    geneExp_fn = os.path.join(RBPexp_dir, sample + '_rpb.csv')
    df = pd.read_csv(geneExp_fn, index_col=0, sep='\t').transpose()
    vec = df.values.flatten() / geneExp_absmax 
    if phase == 'predict':
        vec[vec > 1] = 1
    mat = np.tile(vec, (nrow, 1))

    return mat

def read_splicing_amount(sample1, sample2, splicing_max, eid_list, phase='train', splicing_dir=''):

    splicing_amount_max_eidlist = splicing_max.loc[eid_list]
    
    splicing_amount_fn1 = os.path.join(splicing_dir, sample1+'.output')
    splicing_amount_fn2 = os.path.join(splicing_dir, sample2+'.output')
    
    df1 = pd.read_csv(splicing_amount_fn1, index_col=0, sep='\t')
    df1 = df1.loc[eid_list, 'splicing_amount']
    df1 = df1/splicing_amount_max_eidlist['max_splicing_amount']
    if phase == 'predict':
        df1[df1 > 1] = 1
    df1 = np.asarray(df1)
    df1 = df1.reshape(df1.shape[0], 1)
    
    df2 = pd.read_csv(splicing_amount_fn2, index_col=0, sep='\t')
    df2 = df2.loc[eid_list, 'splicing_amount']
    df2 = df2/splicing_amount_max_eidlist['max_splicing_amount']
    if phase == 'predict':
        df2[df2 > 1] = 1
    df2 = np.asarray(df2)
    df2 = df2.reshape(df2.shape[0], 1)
    
    colnames = ['splicing_amount'] ##

    return df1, df2, colnames

def construct_data_from_label(label_fn, sample, seqFeature_df, geneExp_absmax, geneExp_colnames, RBPexp_dir, splicing_max='', splicing_dir='', splicing_amount=False, sep='_', CIRIdeepA=False, phase='train'):
    
    if CIRIdeepA:
        if phase == 'train':
            no_diff, a_higher, b_higher = read_label_fn(label_fn, CIRIdeepA=CIRIdeepA)
            #***direction***#
            no_diff_eid = list(seqFeature_df.index.intersection(no_diff.keys()))
            a_higher_eid = list(seqFeature_df.index.intersection(a_higher.keys()))
            b_higher_eid = list(seqFeature_df.index.intersection(b_higher.keys()))
            eid_list = no_diff_eid + a_higher_eid + b_higher_eid + no_diff_eid + a_higher_eid + b_higher_eid
        elif phase == 'predict':
            eid_list = list(seqFeature_df.index.intersection(read_label_fn_predict(label_fn)))
            
        seqFeature_df = seqFeature_df.loc[eid_list]
        
        
        # trans-feature
        sample_a = sample.split('_')[0]
        sample_b = sample.split('_')[1]
    
        #***direction***#
        if phase == 'train':
            geneExp_a = read_geneExp(sample_a, geneExp_absmax, nrow=int(len(eid_list)/2), RBPexp_dir=RBPexp_dir)
            geneExp_b = read_geneExp(sample_b, geneExp_absmax, nrow=int(len(eid_list)/2), RBPexp_dir=RBPexp_dir)
            geneExp_left = np.vstack((geneExp_a, geneExp_b))
            geneExp_right = np.vstack((geneExp_b, geneExp_a))
        elif phase == 'predict':
            geneExp_a = read_geneExp(sample_a, geneExp_absmax, nrow=len(eid_list), phase='predict', RBPexp_dir=RBPexp_dir)
            geneExp_b = read_geneExp(sample_b, geneExp_absmax, nrow=len(eid_list), phase='predict', RBPexp_dir=RBPexp_dir)
            
        if phase == 'train':
            X = np.hstack((seqFeature_df, geneExp_left, geneExp_right))
            Y = np_utils.to_categorical(
                np.asarray([0.]*len(no_diff_eid) + [1.]*len(a_higher_eid) + [2.]*len(b_higher_eid) + 
                [0.]*len(no_diff_eid) + [2.]*len(a_higher_eid) + [1.]*len(b_higher_eid)))
            rownames = np.asarray([sample + '-' + x for x in eid_list])
            geneExp_simple_colnames = [gene for gene in geneExp_colnames]
            colnames = np.hstack([seqFeature_df.columns.values, geneExp_colnames])
        elif phase == 'predict':
            X = np.hstack((seqFeature_df, geneExp_a, geneExp_b))
            rownames = np.asarray([sample + '-' + x for x in eid_list])
            geneExp_simple_colnames = [gene for gene in geneExp_colnames]
            colnames = np.hstack([seqFeature_df.columns.values, geneExp_colnames])
        
        if phase == 'train':
            return X, Y, rownames, colnames
        elif phase == 'predict':
            return X, rownames, colnames
    
    else:
        if phase == 'train':
            pos, neg = read_label_fn(label_fn, CIRIdeepA=CIRIdeepA)
            pos_eid = list(seqFeature_df.index.intersection(pos.keys()))
            neg_eid = list(seqFeature_df.index.intersection(neg.keys()))
            eid_list = pos_eid + neg_eid
        elif phase == 'predict':
            eid_list = list(seqFeature_df.index.intersection(read_label_fn_predict(label_fn)))
        
        # cis-feature
        seqFeature_df = seqFeature_df.loc[eid_list]
        
        # trans-feature
        sample_a = sample.split(sep)[0]
        sample_b = sample.split(sep)[1]
        geneExp_a = read_geneExp(sample_a, geneExp_absmax, nrow=len(eid_list), phase=phase, RBPexp_dir=RBPexp_dir)
        geneExp_b = read_geneExp(sample_b, geneExp_absmax, nrow=len(eid_list), phase=phase, RBPexp_dir=RBPexp_dir)
        
        # splicing amount 
        splicing_amount_a, splicing_amount_b, splicing_amount_colnames = read_splicing_amount(sample_a, sample_b, splicing_max, eid_list, phase=phase, splicing_dir=splicing_dir)
        
        X = np.hstack((seqFeature_df, geneExp_a, geneExp_b, splicing_amount_a, splicing_amount_b))
        if phase == 'train':
            Y = np.asarray([1.]*len(pos_eid) + [0.]*len(neg_eid))
        
        rownames = np.asarray([sample + '-' + x for x in eid_list])
        geneExp_simple_colnames = [sample + '-' + gene for sample in [sample_a, sample_b] for gene in geneExp_colnames]
        colnames = np.hstack([seqFeature_df.columns.values, geneExp_simple_colnames, ['splicing_a', 'splicing_b']])
        
        if phase == 'train':
            return X, Y, rownames, colnames
        elif phase == 'predict':
            return X, rownames, colnames


def write_current_pred(test_data, y_pred, fn, CIRIdeepA, phase='train'):

    if phase == 'predict':
        if CIRIdeepA:
            with open(fn, 'w') as fo:
                for i in range(len(test_data['rownames'])):
                    fo.write('\t'.join([test_data['rownames'][i], str(y_pred[i, 0]), str(y_pred[i, 1]), str(y_pred[i, 2])]) + '\n')
        else:
            with open(fn, 'w') as fo:
                for i in range(len(test_data['rownames'])):
                    fo.write('\t'.join([test_data['rownames'][i], str(y_pred[i, 0])]) + '\n')
    elif phase == 'train':
        if CIRIdeepA:
            with open(fn, 'w') as fo:
                for i in range(len(test_data['rownames'])):
                    fo.write('\t'.join([test_data['rownames'][i], str(test_data['Y'][i, 0]), str(test_data['Y'][i, 1]), str(test_data['Y'][i, 2]), str(y_pred[i, 0]), str(y_pred[i, 1]), str(y_pred[i, 2])]) + '\n')
        else:
            with open(fn, 'w') as fo:
                for i in range(len(test_data['rownames'])):
                    fo.write('\t'.join([test_data['rownames'][i], str(test_data['Y'][i]), str(y_pred[i, 0])]) + '\n')

    return

def plot_auroc(roc_val, roc_test, loss_training, loss_eval, loss_test, test_freq, odir):
    
    # plot
    plt.figure(figsize=(6,6))
    plt.title('Evaluation/Test ROC')
    plt.plot(range(1, len(roc_val) + 1), roc_val)
    plt.plot(range(test_freq, len(roc_val) + 1, test_freq), roc_test)
    plt.ylabel('Auroc')
    plt.xlabel('Step')
    plt.savefig(os.path.join(odir, 'evaluation test roc.png'))
    
    # plot loss
    plt.figure(figsize=(6,6))
    plt.title('Training/Evaluation/Test loss')
    plt.plot(range(1, len(loss_training) + 1), loss_training)
    plt.plot(range(1, len(loss_eval) + 1), loss_eval)
    plt.plot(range(test_freq, len(loss_eval) + 1, test_freq), loss_test)
    plt.ylabel('Loss')
    plt.xlabel('Step')
    plt.savefig(os.path.join(odir, 'loss_training_eval_test.png'))
    return

def read_splicing_amount_max(fn):
    df = pd.read_csv(fn, index_col=0, sep='\t')
    return df

def eval_continuous_label_auc(y_true, y_pred, outputdir, prefix):
    auroc = metrics.roc_auc_score(y_true, y_pred)
    fpr, tpr, threshold = metrics.roc_curve(y_true, y_pred)
    np.save(os.path.join(outputdir, 'y_true.npy'), y_true)
    np.save(os.path.join(outputdir, 'y_pred.npy'), y_pred)
    roc_auc = metrics.auc(fpr, tpr)
    plt.figure(figsize=(6,6))
    plt.title('Test AUROC')
    plt.plot(fpr, tpr, 'b', label = 'Test AUC = %0.3f' % roc_auc)
    plt.legend(loc = 'lower right')
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.savefig(os.path.join(outputdir, prefix+'.png'))
    return

def read_roc(odir):
    fn = os.path.join(odir, 'validation_test_auroc.txt')
    roc_val = []
    roc_test = []
    with open(fn, 'r') as f:
        for line in f:
            if line.startswith('v'):
                continue
            content = line.strip().split('\t')
            roc_val.append(content[0])
            if content[1] != 'none':
                roc_test.append(content[1])
    return roc_val, roc_test

def train_model(fn, odir, geneExp_absmax, seqFeature, RBP_dir, splicing_max='', splicing_dir='', splicing_amount=True, CIRIdeepA=False):
    
    if not os.path.isdir(odir):
        os.mkdir(odir)
    
    print('output directory: %s' % odir)
    print('input train list: %s' % fn)
    train_list = read_train_list(fn)
    
    m = clf(**config.architecture, splicing_amount=splicing_amount, CIRIdeepA=CIRIdeepA)
    test_data = {'X':np.empty((0, m.n_in), dtype='float32'), 'Y':np.asarray([], dtype='float32'), 'rownames':np.asarray([], dtype='str')}
    test_data_lst = {'X':[], 'Y':[], 'rownames':[]}
    
    print('Loading features...')
    seqFeature_df = read_sequence_feature(seqFeature)
    geneExp_absmax, geneExp_colnames = read_geneExp_absmax(geneExp_absmax)
    
    if splicing_amount:
        splicing_max = read_splicing_amount_max(splicing_max)
    else:
        splicing_max = []
    
    print('Training start...')
    data_batch_generator = generate_data_batch(
        train_list, test_data_lst, seqFeature_df, geneExp_absmax, geneExp_colnames, odir, RBP_dir, splicing_dir=splicing_dir, splicing_max=splicing_max, splicing_amount=splicing_amount, CIRIdeepA=CIRIdeepA)
    
    # initiate
    roc_val = []
    roc_test = []
    loss_train_lst = []
    acc_train_lst = []
    loss_eval_lst = []
    acc_eval_lst = []
    loss_test_lst = []
    acc_test_lst = []
    best_roc = 0
    best_loss = float('Inf')
    best_loss_mean = float('Inf')
    patience_loss_mean = 0
    patience_loss_mean_limit = 10
    patience_limit = config.patience_limit
    n_batch = 0
    patience = 0
    step = 0
    test_freq = math.ceil(len(train_list.index)/4)
    print('step', 'auroc_eval', 'aupr_eval', 'loss_eval', 'acc_eval', 'loss_train', 'acc_train', sep='\t', file=open(odir + '/roc_pr_losseval_losstrain.log', 'a+'))
    print('step', 'patience', 'auroc_test', 'aupr_test', 'loss_test', 'acc_test', sep='\t', file=open(odir + '/roc_pr_loss_test.log', 'a+'))
    
    # training start
    for data_batch in data_batch_generator:
        
        step += 1
                    
        X = data_batch['X']
        
        Y = data_batch['Y']
        rownames = data_batch['rownames']
        
        # train and evaluation
        inputdata = {'X':X, 'Y':Y, 'rownames':rownames}
        X_train, X_val, Y_train, Y_val, rownames_train, rownames_val = split_eval_train_data(inputdata, n_fold=config.n_fold)
        
        inputdata = {'X_train':X_train, 'X_val':X_val, 'Y_train':Y_train, 'Y_val':Y_val, 'rownames_train':rownames_train, 'rownames_val':rownames_val}
        print('step %i optimization start.' % step)
        best_metrics = m.fit(inputdata)
        if not CIRIdeepA:
            loss_train_lst.append(best_metrics['loss_train'][0])
            acc_train_lst.append(best_metrics['loss_train'][1])
            loss_eval_lst.append(best_metrics['loss_eval'][0])
            acc_eval_lst.append(best_metrics['loss_eval'][1])
            print(step, best_metrics['roc'], best_metrics['aupr'], best_metrics['loss_eval'][0], best_metrics['loss_eval'][1], best_metrics['loss_train'][0], best_metrics['loss_train'][1], sep='\t', file=open(odir + '/roc_pr_losseval_losstrain.log', 'a+'))
            roc_val.append(best_metrics['roc'])
        else:
            loss_train_lst.append(best_metrics['loss_train'])
            acc_train_lst.append(best_metrics['acc_train'])
            loss_eval_lst.append(best_metrics['loss_val'])
            acc_eval_lst.append(best_metrics['acc_val'])
            print('%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' % 
                (step, best_metrics['precision_score_macro'], best_metrics['precision_score_weighted'], 
                best_metrics['recall_macro'], best_metrics['recall_weighted'], 
                best_metrics['F1_macro'], best_metrics['F1_weighted'], 
                best_metrics['loss_val'], best_metrics['acc_val'], 
                best_metrics['loss_train'], best_metrics['acc_train']), 
                file=open(odir + '/roc_pr_losseval_losstrain.log', 'a+'))

        # metrics of test dataset
        if step % test_freq == 0:
                
            if step <= test_freq:
                test_data_lst['X'] = np.concatenate(test_data_lst['X'])
                test_data_lst['Y'] = np.concatenate(test_data_lst['Y'])
                test_data_lst['rownames'] = np.concatenate(test_data_lst['rownames'])
            
            loss_test = m.model.evaluate(test_data_lst['X'], test_data_lst['Y'], verbose=1)
            Y_pred = m.predict(test_data_lst)
            
            if not CIRIdeepA:
                auroc_test = metrics.roc_auc_score(test_data_lst['Y'], Y_pred)
                aupr_test = metrics.average_precision_score(test_data_lst['Y'], Y_pred)
                roc_test.append(auroc_test)
                loss_test_lst.append(loss_test[0])
                acc_test_lst.append(loss_test[1])
                    
                # write predict result
                write_current_pred(test_data_lst, Y_pred, os.path.join(odir, 'step_{0}_pred.txt'.format(step)), CIRIdeepA, phase='train')
                
                # plot auroc
                plot_auroc(roc_val, roc_test, loss_train_lst, loss_eval_lst, loss_test_lst, test_freq, odir)
            
            else:
                Y_predict_categorical = (Y_pred == Y_pred.max(axis=1, keepdims=True)).astype('float')
                precision_score_macro = metrics.precision_score(test_data_lst['Y'], Y_predict_categorical, average="macro")
                precision_score_weighted = metrics.precision_score(test_data_lst['Y'], Y_predict_categorical, average="weighted")
                recall_macro = metrics.recall_score(test_data_lst['Y'], Y_predict_categorical, average='macro')
                recall_weighted = metrics.recall_score(test_data_lst['Y'], Y_predict_categorical, average="weighted")
                F1_macro = metrics.f1_score(test_data_lst['Y'], Y_predict_categorical, average='macro')
                F1_weighted = metrics.f1_score(test_data_lst['Y'], Y_predict_categorical, average='weighted')
                loss_test_lst.append(loss_test[0])
                acc_test_lst.append(loss_test[1])
                
                write_current_pred(test_data_lst, Y_pred, os.path.join(odir, 'step_{0}_pred.txt'.format(step)), CIRIdeepA, phase='train')
                
                # plot auroc
                plot_auroc(acc_eval_lst, acc_test_lst, loss_train_lst, loss_eval_lst, loss_test_lst, test_freq, odir)
            
            # compare the metrics
            if loss_test[0] < best_loss:
                best_loss = loss_test[0]
                m.model.save(os.path.join(odir, str(step)+'_model_weight.h5'))
                patience = 0
            else:
                patience += 1
                m.model.save(os.path.join(odir, str(step)+'_model_weight.h5'))
                if patience > patience_limit:
                    print('patience > patience_limit. Early Stop.')
                    break
            del m
            K.clear_session()
            m = clf(**config.architecture, splicing_amount=splicing_amount, CIRIdeepA=CIRIdeepA)
            m.model.load_weights(os.path.join(odir, str(step)+'_model_weight.h5'))
            
            # print test result 
            print('*' * 50)
            print('test result:')
            if not CIRIdeepA:
                print('step %i: auroc-%.3f; aupr-%.3f; best_loss=%.3f    patience-%i' % (step, auroc_test, aupr_test, best_loss, patience))
                print('*' * 50)
                print(step, patience, auroc_test, aupr_test, loss_test[0], loss_test[1], sep='\t', file=open(odir + '/roc_pr_loss_test.log', 'a+'))
            else:
                print('step %i: acc-%.3f; loss=%.3f; best loss=%.3f; patience-%i' % (step, loss_test[1], loss_test[0], best_loss, patience))
                print('*' * 50)
                print('%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' % 
                    (step, patience, precision_score_macro, precision_score_weighted, 
                    recall_macro, recall_weighted, 
                    F1_macro, F1_weighted, 
                    loss_test[0], loss_test[1]), 
                    file=open(odir + '/roc_pr_loss_test.log', 'a+'))
    
    m.model.save(os.path.join(odir, 'final_model_weight.h5'))
    Y_pred = m.predict(test_data)
    auroc_test = metrics.roc_auc_score(test_data['Y'], Y_pred)
    print('auroc after 10 epoch:', auroc_test)
    
    eval_continuous_label_auc(test_data['Y'], Y_pred, odir, 'auroc')
    
    return


def main():
    
    prog = 'CIRI-deep'
    description = 'CIRI-deep predicts differentially spliced circRNAs between biological samples. This program can be used for training and predicting.'
    epilog = 'Zihan Zhou <zhouzihan2018m@big.ac.cn>'
    argparser = ap.ArgumentParser(prog=prog, description=description, epilog=epilog)
    subparser = argparser.add_subparsers(dest='subcommand')
    
    train = subparser.add_parser('train', help='Train CIRI-deep or CIRI-deep(A) using customized data.')
    predict = subparser.add_parser('predict', help='Predict differentially spliced circRNA using CIRI-deep.')
    
    predict.add_argument('-predict_list', help='file to predict.', dest='predict_list', required=True)
    predict.add_argument('-model_path', help='model path.', dest='model_path', required=True)
    predict.add_argument('-outdir', help='output directory.', dest='outdir', required=True)
    predict.add_argument('-geneExp_absmax', help='gene expression max value.', dest='geneExp_absmax', required=True)
    predict.add_argument('-seqFeature', help='cis feature for all circular RNA.', dest='seqFeature', required=True)
    predict.add_argument('-RBP_dir', help='RBP directory.', dest='RBP_dir', required=True)
    predict.add_argument('-splicing_max', help='splicing amount max value.', default='', dest='splicing_max')
    predict.add_argument('-splicing_dir', help='splicing_amount directory.', default='', dest='splicing_dir')
    predict.add_argument('-sep', help='seperate two samples.', default='_', dest='sep')
    predict.add_argument('--CIRIdeepA', action='store_true', dest='CIRIdeepA')

    train.add_argument('-train_list', help='file to predict.', dest='train_list', required=True)
    train.add_argument('-outdir', help='output directory.', dest='outdir', required=True)
    train.add_argument('-geneExp_absmax', help='gene expression max value.', dest='geneExp_absmax', required=True)
    train.add_argument('-seqFeature', help='cis feature for all circular RNA.', dest='seqFeature', required=True)
    train.add_argument('-splicing_max', help='splicing amount max value.', dest='splicing_max', default='')
    train.add_argument('-RBP_dir', help='RBP directory.', dest='RBP_dir', required=True)
    train.add_argument('-splicing_dir', help='splicing_amount directory.', dest='splicing_dir', default='')
    train.add_argument('--CIRIdeepA', action='store_true', dest='CIRIdeepA')
        
    args = argparser.parse_args()
    subcommand = args.subcommand
    
    if subcommand == 'train':
        train_list = args.train_list
        outdir = args.outdir
        geneExp_absmax = args.geneExp_absmax
        seqFeature = args.seqFeature
        splicing_max = args.splicing_max
        splicing_dir = args.splicing_dir
        RBP_dir = args.RBP_dir
        CIRIdeepA = args.CIRIdeepA
        if CIRIdeepA:
            splicing_amount = False
        else:
            splicing_amount = True
            
        train_model(train_list, outdir, geneExp_absmax, seqFeature, RBP_dir, splicing_max=splicing_max, splicing_dir=splicing_dir, splicing_amount=splicing_amount, CIRIdeepA=CIRIdeepA)
    
    if subcommand == 'predict':
        geneExp_absmax = args.geneExp_absmax
        seqFeature = args.seqFeature
        splicing_max = args.splicing_max
        test_fn = args.predict_list
        model_path = args.model_path
        out_dir = args.outdir
        RBP_dir = args.RBP_dir
        splicing_dir = args.splicing_dir
        sep = args.sep
        CIRIdeepA = args.CIRIdeepA
        if CIRIdeepA:
            splicing_amount = False
        else:
            splicing_amount = True
            
        phase = 'predict'
        
        print('Loading file from %s' % seqFeature)
        seqFeature_df = read_sequence_feature(seqFeature)
        
        print('Loading gene expression from %s' % geneExp_absmax)
        geneExp_absmax, geneExp_colnames = read_geneExp_absmax(geneExp_absmax)
        
        if splicing_amount:
            print('Loading splicing amount feature from %s' % splicing_max)
            splicing_max = read_splicing_amount_max(splicing_max)
        
        print('Loading model weight from %s' % model_path)
        clf = load_model(model_path)
        test_list = read_train_list(test_fn)
        
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        
        print('Prediction...')
        for i in range(len(test_list)):
            inputfile = test_list.label_fn[i]
            sample = test_list.index[i]
            X_tmp, rownames_tmp, colnames = construct_data_from_label(inputfile, 
              sample, seqFeature_df, geneExp_absmax, geneExp_colnames, RBP_dir, splicing_max=splicing_max, splicing_amount=splicing_amount, phase=phase, splicing_dir=splicing_dir, sep=sep, CIRIdeepA=CIRIdeepA)
            test_data = {'X':X_tmp, 'rownames':rownames_tmp}
            y_pred = clf.predict([test_data['X']])
            write_current_pred(test_data, y_pred, '%s/%s.txt' % (out_dir, sample), CIRIdeepA, phase='predict')
            
    return 

if __name__ == '__main__':
    main()