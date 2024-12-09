# $(\ell,\ell)$-isogeny from a product of ellptic curves

このパッケージは[SCIS2025](https://www.iwsec.org/scis/2025/index.html)で吉住(九州大学)が発表したアルゴリズムのSageMathによる実装です. アルゴリズムの詳しい内容は(参加者のみが見れる)予稿を参照します. 

This package contains the implementation of the algorithm in SageMath presented at [SCIS2025](https://www.iwsec.org/scis/2025/index_en.html) by Ryo Yoshizumi. The detailed contents of the algorithm can be found in the proceedings, which are accessible only to participants.


## Content 

In this package, there is functions to compute $(\ell,\ell)$-isogeny from a product of ellptic curves. The main module is `func_isogeny.py`. For the overall flow, please refer to `test.py`.

This code is written in [SageMath](https://www.sagemath.org).



## One Example


By writing the following command, you can compute one example of $(\ell,\ell)$-isogeny from $E_0\times E_0$. 

``` 
% sage example.py
 ```

Here, the base field is $\mathbb{F}_{p^2}$ where $p=276154505650672190920223$, $\ell=11$, and $E_0 : y^2=x^3+x$.

The kernel $K\subset (E_0\times E_0)[\ell]$  is ***randomly*** generalted, and we write  a basis of $K$ by $e_1=(e_1^{(1)},e_1^{(2)}),e_2=(e_2^{(1)},e_2^{(2)})$.

Evaluating point is $x=(x^{(1)},0_{E_0})$ where $x^{(1)}\in E_0$ is a random point.

The output is as follows:
 ```
 p: 276154505650672190920223

 ell: 11

 E_0:
 Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 276154505650672190920223^2

 e1=(e1^(1),e1^(2)):
 ((210831735326664260671440*z2 + 87198234199232688535215 : 133255005542521365746194*z2 + 164977575746837862238202 : 1), (83375910444643963132697*z2 + 200948769993679858502156 : 144292860645864077727680*z2 + 148419487573216577660120 : 1))

 e2=(e2^(1),e2^(2)):
 ((98830049968673606594523 : 242317345512027435674284*z2 + 49548615497841843339319 : 1), (177324455681998584325700 : 216601895301696610344217 : 1))

 x=(x^(1),0):
 ((106189123109636747922251*z2 + 248968408709573238167579 : 151118498536791259433057*z2 + 87307954336256715411728 : 1), (0 : 1 : 0))

 theta-null point of the codomain:
 [112027156618648925977355*z2 + 147438371256409604632158, 69068842046189523671338*z2 + 222257352160280440335182, 275376810320619078414513*z2 + 183273802569472479972437, 86875977424383747210632*z2 + 210469468339113327744966, 1]

 theta coordinate of the image f(x):
 [51118548117130760667387*z2 + 22524482140569545457684, 179878823554441321130688*z2 + 205994074250712649272564, 198607624210760764176903*z2 + 126880054694919192506851, 103634261912850587289929*z2 + 224536590004812764733705, 1]
 ```



## Author
- Name: Ryo Yoshizumi 
- Affiliation: Kyushu University
- Email: yoshizumi.ryo.483@s.kyushu-u.ac.jp





