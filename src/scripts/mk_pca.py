import pandas as pd
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt
%matplotlib inline

rnaseq_samples_mat = DATA + 'interim/TPM.txt'
df = pd.read_csv(rnaseq_samples, sep='\t')

crit = df.sum(axis = 1) > 
dat = df[crit] # selects True rows

X = dat[features].values
X = StandardScaler().fit_transform(X)
colors = {'ill':'red',
          'both':'green',
          'cgi':'blue',
          'none':'black'}
status_colors = [colors[x] for x in dat['kaviar_status']]
print('here')
pca = PCA(n_components=2)
p = pca.fit(X)
print(p.explained_variance_ratio_)
X_r = p.transform(X)

plt.figure()
#for x_y, status in zip(X_r, status):
plt.scatter(X_r[:, 0], X_r[:, 1], alpha=.4, c=status_colors)

##########
pca = PCA()
X_r = pca.fit(X).transform(X)

print(pca.explained_variance_ratio_)

plt.figure()
colors = ['navy', 'darkorange']
lw = 2

for color, i, target_name in zip(colors, [0, 1], ['cgi', 'both']):
    plt.scatter(X_r[y == i, 0], X_r[y == i, 1], alpha=.5, color=color,
                label=target_name)
plt.legend(loc='best', shadow=False, scatterpoints=1)
plt.title('PCA of RNA seq samples')

plt.show()
