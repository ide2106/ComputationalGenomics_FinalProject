import pyreadr
import math
from scipy import stats
import csv
import matplotlib.pyplot as plt
import numpy as np


class MR_pipeline():

    def __init__(self, gene_f, exp_f, outcome_f, pval=0.01, print_results=True, print_graphs=False, check_missing=None, exposure_col=0):
        self.gene_file = gene_f
        self.exposure_file = exp_f
        self.outcome_file = outcome_f
        self.giv_vs_exp_pval = pval
        self.verbose = print_results
        self.graphs = print_graphs
        self.check_missing = check_missing
        self.exp_col = exposure_col


    def full_analysis(self):

        try:
            if self.exposure_file.endswith('.csv'):
                exp_dicts = self.csv_data(self.exposure_file, is_float=True)
            elif self.exposure_file.endswith('.RData'):
                exp_dicts = self.rdata_data(self.exposure_file, is_float=True)
            else:
                print('Error: cannot read %s, must input in .csv or .RData files' % self.exposure_file)
                return 0
        except:
            print('Error: %s not formatted correctly' % self.exposure_file)
            return 0

        exp = exp_dicts[self.exp_col] #takes first column in exposure file - SHOULD ONLY HAVE 1

        try:
            if self.gene_file.endswith('.csv'):
                giv_dicts = self.csv_data(self.gene_file, is_int=True)
            elif self.gene_file.endswith('.RData'):
                giv_dicts = self.rdata_data(self.gene_file, is_int=True)
            else:
                print('Error: cannot read %s, must input in .csv or .RData files' % self.gene_file)
                return 0
        except:
            print('Error: %s not formatted correctly' % self.gene_file)
            return 0

        max_corr = 0
        best_giv = None
        giv_vs_exp = []
        for i in range(0,len(giv_dicts)):
            giv = giv_dicts[i]
            slope, intercept, R, pr, stderr = self.snp_vs_exp(giv, exp)[1] # [1] = regression results (0 = chi)
            if pr < self.giv_vs_exp_pval and abs(slope) > max_corr:
                max_corr = abs(slope)
                best_giv = i
                giv_vs_exp = [slope, intercept, R, pr, stderr]

        if best_giv is not None:
            if self.verbose:
                print('The %s gene in %s was the strongest predictor of exposure, with regression slope '
                      '%f and intercept %f, p val:' % (best_giv, self.gene_file, giv_vs_exp[0], giv_vs_exp[1]), giv_vs_exp[3])
                print()
        else:
            print('Error: No gene in %s predicted given exposure')
            return 0
        try:
            if self.outcome_file.endswith('.csv'):
                outcome_data = self.csv_data(self.outcome_file, is_float=True)
            elif self.outcome_file.endswith('.RData'):
                outcome_data = self.rdata_data(self.outcome_file, is_float=True)
            else:
                print('Error: cannot read %s, must input in .csv or .RData files' % self.outcome_file)
                return 0
        except:
            print('Error: %s not formatted correctly' % self.outcome_file)
            return 0

        if self.graphs:
            self.SNP_vs(exp, giv, 'Exposure by GIV')

        good = []
        for i in range(0, len(outcome_data)):
            data = self.given_exp(giv, exp, outcome_data[i], np.unique(list(exp.values())), np.unique(list(exp.keys())))

            none_count = 0
            for k in outcome_data[i].keys():
                if outcome_data[i][k] is None:
                    none_count += 1

            independent_flag = True
            for d in data:
                if d[3] < 0.01:
                    independent_flag = False

            if independent_flag:
                if self.check_missing is None or none_count < self.check_missing:
                    good.append(i)
                    if self.verbose:
                        print('MR assumptions hold on %d column in %s' % (i, self.outcome_file))
                        print('given exp:')
                        print(data)
                        print()

                else:
                    if self.verbose:
                        print('MR assumptions DO NOT hold on %d column in %s' % (i, self.outcome_file))
                        print('None count: %d' %none_count)
                        print('given exp:')
                        print(data)
                        print()

        if self.verbose:
            print('Outcome data to use for MR')
            for i in good:
                print(i)

        for i in good:
            slope, intercept, R, pr, stderr = self.MR(giv, exp, outcome_data[i])
            if self.verbose:
                print('MR results for', i)
                if pr < 0.01:
                    print('Causal relationship suggested')
                    print('Estimated exposure predicts outcome with slope: %f, intercept: %f, p-val: ' %(slope, intercept), pr)
                    print()
                else:
                    print('Cannot predict causal relationship with this data and p-val < 0.01')
                    print('Estimated exposure predicts outcome with slope: %f, intercept: %f, p-val: ' %(slope, intercept), pr)
                    print()


    def csv_data(self, filename, headers=False, is_int=False, is_float=False):

        to_return = []
        with open(filename, encoding='utf-8-sig') as f:
            data = csv.reader(f)
            count = 0
            for row in data:
                if count == 0:
                    count += 1
                    for x in range(0, len(row)-1):
                        D = {}
                        to_return.append(D)
                    if headers:
                        continue

                for i in range(0, len(row)):
                    if i == 0:
                        continue
                    else:
                        if is_int or is_float:
                            if row[i] != 'NA':
                                if is_int:
                                    to_return[i - 1][int(row[0])] = int(row[i])
                                else:
                                    to_return[i - 1][int(row[0])] = float(row[i])
                            else:
                                to_return[i - 1][int(row[0])] = None
                        else:
                            to_return[i-1][int(row[0])] = row[i]
        return to_return


    def rdata_data(self, filename, is_int=False, is_float=False):

        to_return = []
        data = pyreadr.read_r(filename)

        for k in data.keys():
            for r in range(0, len(data[k].index)):
                row = data[k].iloc[r, :]
                if r == 0:
                    for i in range(0, len(row)-1):
                        D = {}
                        to_return.append(D)

                for j in range(1, len(row)):
                    if is_int or is_float:
                        try:
                            if not math.isnan(float(row[j])):

                                if is_float:
                                    to_return[j - 1][int(row[0])] = float(row[j])
                                else:
                                    to_return[j - 1][int(row[0])] = float(row[j])
                        except:
                            if len(row[j]) != 0:
                                to_return[j - 1][int(row[0])] = float(row[j][1:])
                    else:
                        to_return[j-1][int(row[0])] = row[j]

        return to_return


    def linear_regression(self, SNP_dict, exp_dict, keys, p=False):
        x = []
        y = []
        for k in keys:
            if k in SNP_dict and k in exp_dict and SNP_dict[k] is not None and exp_dict[k] is not None:
                x.append(SNP_dict[k])
                y.append(exp_dict[k])
        if len(x) != 0 and len(y) != 0:
            slope, intercept, R, pr, stderr = stats.linregress(x,y)
        else:
            print('Error: dicts do not share unique IDs for individuals')
            print(SNP_dict)
            print(exp_dict)
            print()
            slope, intercept, R, pr, stderr = [0,0,0,0,0]
        if p:
            print(slope, intercept, R, pr, stderr)

        return slope, intercept, R, pr, stderr


    def MR(self, gene, exposure, outcome):
        exp_bar = {}
        slope, intercept, R, pr, stderr = self.linear_regression(gene, exposure, list(gene.keys()))
        for RID in gene.keys():
            exp_bar[RID] = intercept + slope* gene[RID]

        print(exposure)
        print(exp_bar)

        slope,intercept,R, pr, stderr = self.linear_regression(exp_bar, outcome, list(exp_bar.keys()))
        #MR_graph(exp_bar, outcome)
        return slope, intercept, R, pr, stderr


    def given_exp(self, SNP_dict, exp_dict, outcome_dict, exp_key, keys):
        # first 3 dicts map RID to values
        # exp key maps vals in exp dict to True or False

        SNP_dicts = {}
        outcome_dicts = {}
        for x in exp_key:
            SNP_dicts[x] = {}
            outcome_dicts[x] = {}

        for k in keys:
            if k in SNP_dict and k in outcome_dict:
                SNP_dicts[exp_dict[k]][k] = SNP_dict[k]
                outcome_dicts[exp_dict[k]][k] = outcome_dict[k]

        to_return = []
        for x in exp_key:
            data = self.linear_regression(SNP_dicts[x], outcome_dicts[x], list(SNP_dict.keys()))
            to_return.append(data)

        return to_return


    def snp_vs_exp(self, gene, exposure):

        x_vals = list(gene.values())
        x_vals = [i for i in x_vals if i is not None]
        X = np.unique(x_vals)

        y_vals = list(exposure.values())
        y_vals = [i for i in y_vals if i is not None]
        Y = np.unique(y_vals)

        dist = np.zeros((len(X), len(Y)))
        count= 0
        for RID in gene.keys():
            if RID in exposure and gene[RID] is not None and exposure[RID] is not None:
                dist[list(X).index(gene[RID]), list(Y).index(exposure[RID])] += 1
                count +=1

        expect = np.full(dist.shape, count/dist.size)

        to_ret = []
        to_ret.append(stats.chisquare(dist, expect))
        to_ret.append(self.linear_regression(gene, exposure, gene.keys()))

        return to_ret


    def exposure_vs_outcome(self, adas, trails):

        x = []
        y = []
        points = []

        for k in adas.keys():
            x.append(adas[k])
            y.append(trails[k])
            points.append([adas[k], trails[k]])

        plt.clf()
        plt.scatter(x, y)
        plt.show()


    def SNP_vs(self, adas, apoe, title):
        apoe0 = []
        apoe1 = []
        apoe2 = []

        for k in adas.keys():
            if apoe[k] == 0 and adas[k] is not None:
                apoe0.append(adas[k])

            if apoe[k] == 1 and adas[k] is not None:
                apoe1.append(adas[k])

            if apoe[k] == 2 and adas[k] is not None:
                apoe2.append(adas[k])

        plt.ylabel('Score Distribution')
        plt.xlabel('APoE4 Allele Count')
        plt.title(title)
        plt.boxplot([apoe0, apoe1, apoe2], positions=[0,1,2])
        plt.show()


    def MR_graph(self, exp_est, outcome):

        x = []
        y = []

        for k in exp_est.keys():
            if k in outcome and exp_est[k] is not None and outcome[k] is not None:
                x.append(exp_est[k])
                y.append(outcome[k])

        plt.ylabel('Outcome')
        plt.xlabel('Predicted Exposure')
        plt.title('MR Analysis ')
        plt.plot(x,y, marker='o', linestyle='')
        plt.show()


