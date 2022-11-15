# BSS vs the rest

Some take aways (at least for me) from working on the issues raised by the reviewers:

## I. MCC vs F1: in the range of $0$ and $\alpha * 0.1348845$
In our simulation study the F1 and MCC values are very similar in all settings. Why does this happen? Especially with MCC ranging in theory from -1 (complete disaster with respect to classification) to +1 (overachiever) while F1 has a range from 0 to 1?
I argue that the similarity between F1 and MCC values in our simulation study comes from the excess of negatives which are common in high dimensional variable selection problems. To be more specific, MCC is defined as

$$\frac{(TP * TN) - (FP * FN)}{\sqrt{(TP + FP) * (TP + FN) * (TN + FP) * (TN +FN)}}$$
 	
while F1 is defined as

$$\frac{2TP}{2TP+FP+FN}.$$ 

Even if MCC can become negative in theory (if $(FP * FN)>(TP * TN)$ ) this is very unlikely to happen in a variable selection problem where one assumes *sparsity* ( $s \ll p$ like in our settings) and small selected subsets ( $k$ ). In these settings the *number of negatives dominates the number of positives* strongly. As long as a selected subset size is rather small compared to the number of variables ( $k \ll p$ ) we will always label most of the variables correctly as negative just by chance/simply because we have to. Consequently, even if all selected variables are wrongly labeled as positive, the total number of FP and FN is small compared to TN. Thus, as long as we will label at least one variable correctly as positive ( $TP \geq 1$ ) we will always have a positives numerator. In our simulation in nearly every run at least one variable is correctly identified as a positive, hence, we have only a handful of negative MCC values.

Now, let us rewrite MCC as 
$$\frac{(TP * TN)}{\sqrt{(TP + FP) * (TP + FN) * (TN + FP) * (TN +FN)}} -$$ 
$$\frac{(FP * FN)}{\sqrt{(TP + FP) * (TP + FN) * (TN + FP) * (TN +FN)}}.$$

Assuming sparsity implies having a rather big TN compared to FP, FN and TP for small $k$. 

E.g.: in our medium dimensional setting we have $p=500$ variable of which $s=10$ are positives (true direct predictors) and the subset sizes range over $k=1,...,15$ for BSS (of course for Enet and Lasso we also have bigger subset sizes but for the moment we are only interested in the smaller ones). This means $TP \leq 10$, $FP \leq 15$, $FN \leq 10$ and $TN \geq 475$.

So it might be fair to say that $(TN + FP) \approx (TN +FN) \approx TN$ is not too rough approximation. Hence, the second term of the previous MCC formulation is rather small in settings similar to ours since it is dominated by 
$$\sqrt{(TN + FP) * (TN +FN)} \approx \sqrt{TN} * \sqrt{TN} = TN$$ 
in the denominator. Therefore I will neglect this term for my further argumentation. 

To give an impression of the highest possible value for the second term in our settings I will use the medium setting for an numerical example. We have $s=10$ true direct predictors (positives), a maximum subset size of $k=15$ and $p=500$ variables. If a selected subset in this setting is completely wrong, i.e. it consists only of FPs, we have FP=15 and FN =10 while still having 475 TN. In this case the second term will be not higher than 0.0251233, its highest possible value (in the high dimensional setting with $p=1000$ it is even smaller, i.e. 0.002).

Using $\sqrt{(TN + FP) * (TN +FN)} \approx TN$ also in the first term of the MCC we get 
$$\frac{(TP * TN)}{\sqrt{(TP + FP) * (TP + FN) * (TN + FP) * (TN +FN)}} \approx \frac{(TP * TN)}{\sqrt{(TP + FP) * (TP + FN)} * TN} = \frac{TP}{\sqrt{(TP + FP) * (TP + FN)}}.$$ 
Let's call this MCC approximation $MCC_{oTN}$ since it omits the TN.

Based on this formulation the F1 score and the $MCC_{oTN}$ both consist only of TP, FP and FN. Since a selected subset size is $k=TP+FP$ and the actual number of positives (true direct predictors) is $s=TP+FN$ we can reformulate the F1-score and $MCC_{oTN}$:
$$F_1=\frac{2TP}{2TP+FP+FN}=\frac{2TP}{k+s}=\frac{TP}{0.5*(k+s)}$$
and
$$MCC_{oTN} = \frac{TP}{\sqrt{k*s}}.$$

Since $k,s \geq 1$ we have $0.5 * (k+s) \geq \sqrt{k*s}$ and therefore 
$$MCC_{oTN} \geq F1.$$ 
For the special case $k=s$ we have $MCC_{oTN} = F1$. 

So far we know that the $MCC_{oTN}$ will always be at least as high as the F1 score (assuming an excess of TN). But we still want to know the range of their difference. To answer this question we first have to take the TP in numerators into account. For $k < s$ its highest possible value is $k$ while for $k \geq s$ it is $s$. Therefore we define

$$\Delta(k,s) := \frac{\min(k,s)}{\sqrt{k+s}} - \frac{\min(k,s)}{0.5*(k+s)}$$

as the difference with the highest number of possible TP for a given subset size $k$. NOTE: we allow also the case $k>s$.

While it might be possible to find the range of $\Delta(k,s)$ algebraically I used Wolfram Alpha to find it (numerically) to speed things up:  
$$0 \leq \Delta(k,s) \leq 0.1348845$$
We now have fixed upper and lower limits for the difference between the two performance measures. And it is not really big. And it gets even smaller since we have only looked at the case when all selected variables are TPs ( $k < s$ ) or all positives have been selected ( $ k \geq s$ ), i.e. the highest possible numerator. As soon as the number of TPs gets smaller the upper limit for the range gets smaller too.

Augmenting $\Delta(k,s)$ with $0 < \alpha \leq 1$ to denote the proportion of possible TPs found for a given subset size $k$ we get
$$\Delta(k,s, \alpha) := \frac{\alpha \min(k,s)}{\sqrt{k+s}} - \frac{\alpha \min(k,s)}{0.5*(k+s)}.$$
From this it follows directly 
$$0 \leq \Delta(k,s, \alpha) \leq \alpha * 0.1348845.$$



