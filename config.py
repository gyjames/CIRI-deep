import os

test_prop = 0.01
random_seed_test = 12345
random_seed_eval = 12345
n_fold = 10
pos_prop = 0.7
class_weight = {0:1, 1:1.2, 2:1.2}

architecture = {
    'n_in': 6523, 
    'hidden': [1200, 500, 300, 200], 
    'drop_out': [0, 0.5, 0.3, 0.2, 0.1], 
    'batch_size': 512, 
    'learning_rate': 0.001
    }

test_freq = 6000
patience_limit = 200

n_epoch = 100

