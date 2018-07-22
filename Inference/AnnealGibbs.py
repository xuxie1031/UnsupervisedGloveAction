'''
Annealing process of fine tuning the learned clustering labels
via Gibbs Sampler
'''

import nltk
import scipy.io as sio
from scipy.stats import multivariate_normal
import numpy as np
import math
import random
import os
import argparse


# load grammar
def read_induced_grammar(path):
    with open(path) as f:
        rules = [rule.strip() for rule in f.readlines()]
        grammar = nltk.PCFG.fromstring(rules)
    	return grammar


# read and parse inputs
def input_data(data_name, label_name, gaussian_name):
	data=sio.loadmat(data_name)
	fvec=data['wrist_vec']

	data=sio.loadmat(label_name)
	labels=data['Cids']
	labels=labels[0,:]

	gaussians=sio.loadmat(gaussian_name)
	label_mean=gaussians['means']
	label_cov=gaussians['covs']

	return fvec, labels, label_mean, label_cov	


# calculate loglikelihood table
def calc_loglikelihood_table(data, num_l, label_mean, label_cov):
	small_prob=1e-200

	num=np.shape(data)[0]
	dim=np.shape(data)[1]
	loglikelihood_table=np.zeros((num, num_l), dtype=np.float64)
	
	for i in range(num):
		for l in range(num_l):
			mean=label_mean[l,:]
			cov=label_cov[l,:].reshape(dim,dim).T
			prob=multivariate_normal.pdf(data[i,:], mean=mean, cov=cov)
			if prob<small_prob:
				prob=small_prob
			loglikelihood_table[i,l]=math.log(prob)

	np.savetxt('loglikelihoods.csv', loglikelihood_table, delimiter=",")	
	return loglikelihood_table


# calculate prior p(pg)
def calc_logprior(grammar, tokens):
	invalid_prob=1e-20

	viterbi_parser=nltk.ViterbiParser(grammar)

	try:
		v_parses=viterbi_parser.parse_all(tokens)
		if v_parses:
			prob=reduce(lambda a, b: a+b.prob(), v_parses, 0)/len(v_parses)
		else:
			print 'fail parse'
			return math.log(invalid_prob)
	except ValueError:
		return math.log(invalid_prob)

	return math.log(prob)


# calculate likelihood p(Gamma|pg)
def calc_loglikelihood(data, labels, likelihood_table):
	small_prob=1e-200

	log_likelihood=0.0
	n=np.shape(data)[0]
	d=np.shape(data)[1]
	for i in range(n):
		l=labels[i]
		log_likelihood+=likelihood_table[i,l]

	return log_likelihood/1e5


# annealing process
def Gibbs_infer_process(data, labels, loglikelihood_table, grammar, num_l):

	# get segment index position
	def get_segment_pos(labels):
		num_seg=1
		c=labels[0]
		pos=0

		while pos<np.shape(labels)[0]:
			if labels[pos]!=c:
				num_seg+=1
				c=labels[pos]
			pos+=1

		seg_poses=np.zeros((num_seg,2), dtype=np.int32)
		seg_count=0
		start=0
		c=labels[0]
		end=0

		while end<np.shape(labels)[0]:
			if labels[end]!=c:
				seg_poses[seg_count,0]=start
				seg_poses[seg_count,1]=end
				start=end
				c=labels[start]
				seg_count+=1
			end+=1
		seg_poses[seg_count,0]=start
	        seg_poses[seg_count,1]=end
	        seg_count+=1
		return seg_poses

	# get token sequence
	def get_token(labels):
		tokens=[]
		c=labels[0]
		tokens.append(str(c))
		pos=1
		while pos<np.shape(labels)[0]:
			if labels[pos]!=c:
				c=labels[pos]
				tokens.append(str(c))
			pos+=1
		return tokens

	seg_poses=get_segment_pos(labels)
	num_seg=np.shape(seg_poses)[0]

	T=1
	alpha=0.99
	num_iter=100

	for i in range(num_iter):
		for p in range(num_seg):
			start=seg_poses[p,0]
			end=seg_poses[p,1]

			prob_l=np.zeros(num_l, dtype=np.float64)
			for l in range(num_l):
				trial=np.copy(labels)
				np.put(trial, np.arange(start,end), [l])
				tokens=get_token(trial)
				log_prior=calc_logprior(grammar, tokens)
                log_likelihood=calc_loglikelihood(data, trial, loglikelihood_table)
                log_posterior=log_prior+log_likelihood
                print 'log_posterior: propose label{}, log prior{}, log likelihood{}'.format(l, log_prior, log_likelihood)
                prob_l[l]=math.exp(1.0/T*log_posterior)

			cum_prob_l=np.cumsum(prob_l)
            r=random.random()
            cum_prob_l/=cum_prob_l[-1]
			print cum_prob_l
            argl=np.where(cum_prob_l>r)[0][0]
			print 'iter:{}, seg[{},{}], label:{}'.format(i, start, end-1, argl)
			
            np.put(labels, np.arange(start,end), [argl])
			print 'tokens:{}'.format(get_token(labels))
		# break

		T*=alpha

	print 'annealing done'
	return labels	

def main():

	# nlabel = 6
	# data_name = 'sample_hand_data'

	# parse necessary arguments
	argparser = argparse.ArgumentParser(description='AnnelGibbs')
	argparser.add_argument('--nlabel', type=int, default=6, metavar='N',
							help='number of labels before inference')
	argparser.add_argument('--data-name', default='sample_hand_data', metavar='Name',
							help='associated data name')
	args = argparser.parse_args()

	# feature directory
	feature_dir = '../Clustering/hand_feature'
	feature_name = os.path.join(feature_prefix, data_name+'.mat')

	# label directory
	label_dir = '../ACA/ACAbin'
	label_name = os.path.join(label_dir, 'ACA_Cids_C_'+data_name+'.mat')

	# gaussian directory
	gaussian_dir = 'Gaussians/gaussian_params'
	gaussian_name = os.path.join(gaussian_dir, data_name+'_gaussian_params.mat')

	# grammar directory
	grammar_dir = 'induced_grammar'
	grammar_name = os.path.join(grammar_dir, 'parser_'+data_name+'.txt')


	random.seed(0)

	print 'read inputs ...'
	grammar = read_induced_grammar(grammar_name)
	data, labels, label_mean, label_cov=input_data(data_name, label_name, gaussian_name)

	print 'offline calculating loglikelihood terms ...'
	if os.path.exists('loglikelihoods.csv'):
		loglikelihood_table = np.genfromtxt('loglikelihoods.csv', delimiter=',')
	else:
		loglikelihood_table = calc_loglikelihood_table(data, nlabel, label_mean, label_cov)

	print 'Gibbs sampling inference on motion labels ...'
	infer_labels = Gibbs_infer_process(data, labels, loglikelihood_table, grammar, nlabel)

	print 'saving inference results ...'
	save_dir = 'motion_labels'
	if not os.path.exists(save_dir):
		os.mkedirs(save_dir)
	np.savetxt(os.path.join(save_dir, data_name+'.txt'), infer_labels, fmt='%d')

	# input_dir='/home/xuxie/matlab/workplace/glove_action_corl/'
	# bottle_name='12_bottle64_open_palm_tf_convert_merged_successes_proc_all_correct'
	# task_name=bottle_name+'_1'
	# data_name='data/'+task_name+'_fvec.mat'
	# label_name='cluster_result/new_Cids_C_6_'+task_name+'_fvec.mat'
	# gaussian_name='gaussians/'+bottle_name+'.mat'
	# grammar_path='induced_grammar/parser_bottle64.txt'
	# num_l=6

	# random.seed(0)

	# print 'data inputs...'
	# grammar=read_induced_grammar(grammar_path)
	# data, labels, label_mean, label_cov=input_data(input_dir, data_name, label_name, gaussian_name)

	# print 'calculaing loglikelihoods...'
	# if os.path.exists('loglikelihoods.csv'):
	# 	loglikelihood_table=np.genfromtxt('loglikelihoods.csv', delimiter=',')
	# else:
	# 	loglikelihood_table=calc_loglikelihood_table(data, num_l, label_mean, label_cov)

	# print 'annealing...'
	# labels=Gibbs_infer_process(data, labels, loglikelihood_table, grammar, num_l)

	# print 'end. saving...'
	# count=0
	# save_name='label_output'
	# fid=open(save_name,'w')
	# for i in labels.tolist():
	# 	fid.write("%d, " % i)
	# 	#count+=1
	# 	#if count%30==0:
	# 	#	fid.write("\n")
	# 	#	count=0
	# fid.close()

if __name__=='__main__':
	main()
