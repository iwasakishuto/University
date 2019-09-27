var tipuesearch = {"pages":[{"title":"HMMのスケーリング","text":"HMMを実装する際に課題となるのが、forward-algorithm で \\(\\alpha\\) を再帰的に求める際に、 \\(\\alpha\\left(\\mathbf{z}_{n-1}\\right)\\) に \\(p(\\mathbf{z}_n|\\mathbf{z}_{n-1})\\) と \\(p(\\mathbf{x}_n|\\mathbf{z}_n)\\) をかけるため、値が非常に小さくなり、計算機の ダイナミックレンジ を超えてしまうことです。 そこで、ここでは \\(\\alpha\\left(\\mathbf{z}_{n}\\right)\\) と \\(\\beta\\left(\\mathbf{z}_{n}\\right)\\) にスケーリングを施し、それらの値が \\(1\\) のオーダーに止まるようにする手法を説明します。 forward-backward \\(\\alpha,\\beta\\) forward-backward algorithm で用いられていた \\(\\alpha,\\beta\\) は以下のように定義されていました。 $$ \\begin{aligned} \\alpha\\left(\\mathbf{z}_{n}\\right) & \\equiv p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n}, \\mathbf{z}_{n}\\right) & (13.34)\\\\ \\beta\\left(\\mathbf{z}_{n}\\right) & \\equiv p\\left(\\mathbf{x}_{n+1}, \\ldots, \\mathbf{x}_{N} | \\mathbf{z}_{n}\\right) & (13.35) \\end{aligned} $$ Scaling factors まず、スケーリングされた \\(\\alpha,\\beta\\) は以下のように表されます。スケーリングによって、 \\(\\alpha\\) は高々 \\(K\\) 個の変数上の確率分布 、 \\(\\beta\\) は2つの条件付き確率の比 になることがわかります。 $$ \\begin{aligned} \\widehat{\\alpha}\\left(\\mathbf{z}_{n}\\right) &=p\\left(\\mathbf{z}_{n} | \\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n}\\right)=\\frac{\\alpha\\left(\\mathbf{z}_{n}\\right)}{p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n}\\right)} & (13.55)\\\\ \\widehat{\\beta}\\left(\\mathbf{z}_{n}\\right) &=\\frac{p\\left(\\mathbf{x}_{n+1}, \\ldots, \\mathbf{x}_{N} | \\mathbf{z}_{n}\\right)}{p\\left(\\mathbf{x}_{n+1}, \\ldots, \\mathbf{x}_{N} | \\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n}\\right)}=\\frac{\\beta\\left(\\mathbf{z}_{n}\\right)}{p\\left(\\mathbf{x}_{n+1}, \\ldots, \\mathbf{x}_{N} | \\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n}\\right)} & (13.61) \\end{aligned} $$ ここで、これらと \\(\\alpha,\\beta\\) を関連付けるためのスケーリング係数 \\(c\\) を導入します。 $$ c_{n}=p\\left(\\mathbf{x}_{n} | \\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n-1}\\right)\\qquad (13.56) $$ すると、 $$ \\begin{aligned} p\\left(\\mathbf{x}_1,\\ldots,\\mathbf{x}_n\\right) &= p\\left(\\mathbf{x}_n | \\mathbf{x}_1,\\ldots,\\mathbf{x}_{n-1}\\right)\\cdots p\\left(\\mathbf{x}_3 | \\mathbf{x}_1,\\mathbf{x}_{2}\\right)p\\left(\\mathbf{x}_2 | \\mathbf{x}_1\\right)p(\\mathbf{x}_1)\\\\ &= c_n\\cdots c_3c_2c_1\\\\ &= \\prod_{m=1}&#94;{n} c_{m} \\end{aligned}\\qquad (13.58) $$ と展開することができるので、 $$ \\begin{aligned} \\alpha\\left(\\mathbf{z}_{n}\\right)&=p\\left(\\mathbf{z}_{n} | \\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n}\\right) p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n}\\right)=\\left(\\prod_{m=1}&#94;{n} c_{m}\\right) \\widehat{\\alpha}\\left(\\mathbf{z}_{n}\\right) & (13.58)\\\\ \\beta\\left(\\mathbf{z}_{n}\\right)&=p\\left(\\mathbf{x}_{n+1}, \\ldots, \\mathbf{x}_{N} | \\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n}\\right)\\widehat{\\beta}\\left(\\mathbf{z}_{n}\\right) = \\left(\\prod_{m=n+1}&#94;{N} c_{m}\\right) \\widehat{\\beta}\\left(\\mathbf{z}_{n}\\right) & (13.60) \\end{aligned} $$ と対応関係がわかります。 \\(\\gamma,\\xi\\) 続いて、 \\(\\gamma,\\xi\\) と \\(\\widehat{\\alpha},\\widehat{\\beta}\\) の対応関係を求めます。 \\(\\alpha,\\beta\\) \\(\\alpha,\\beta\\) を用いると、以下のように表されていました。 $$ \\begin{aligned} \\gamma\\left(\\mathbf{z}_{n}\\right) &= \\frac{\\alpha\\left(\\mathbf{z}_{n}\\right) \\beta\\left(\\mathbf{z}_{n}\\right)}{p(\\mathbf{X})} & (13.33)\\\\ \\xi\\left(\\mathbf{z}_{n-1}, \\mathbf{z}_{n}\\right) &=\\frac{\\alpha\\left(\\mathbf{z}_{n-1}\\right) p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) p\\left(\\mathbf{z}_{n} | \\mathbf{z}_{n-1}\\right) \\beta\\left(\\mathbf{z}_{n}\\right)}{p(\\mathbf{X})} &(13.43) \\end{aligned}$$ \\(\\widehat{\\alpha},\\widehat{\\beta}\\) 先の対応関係を用いれば、 \\(\\widehat{\\alpha},\\widehat{\\beta}\\) を用いると、 $$ \\begin{aligned} \\gamma\\left(\\mathbf{z}_{n}\\right) &=\\widehat{\\alpha}\\left(\\mathbf{z}_{n}\\right) \\widehat{\\beta}\\left(\\mathbf{z}_{n}\\right) & (13.64)\\\\ \\xi\\left(\\mathbf{z}_{n-1}, \\mathbf{z}_{n}\\right) &=\\left(c_{n}\\right)&#94;{-1} \\widehat{\\alpha}\\left(\\mathbf{z}_{n-1}\\right) p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) p\\left(\\mathbf{z}_{n} | \\mathbf{z}_{-1}\\right) \\widehat{\\beta}\\left(\\mathbf{z}_{n}\\right) & (13.65) \\end{aligned} $$ と表されることがわかります。 Recursion 最後に、再帰式の対応関係も求めます。 \\(\\alpha,\\beta\\) $$ \\begin{aligned} \\alpha\\left(\\mathbf{z}_{n}\\right) &=p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) \\sum_{\\mathbf{z}_{n-1}} \\alpha\\left(\\mathbf{z}_{n-1}\\right) p\\left(\\mathbf{z}_{n} | \\mathbf{z}_{n-1}\\right) & (13.36)\\\\ \\alpha\\left(\\mathbf{z}_{1}\\right) &=p\\left(\\mathbf{x}_{1}, \\mathbf{z}_{1}\\right)=p\\left(\\mathbf{z}_{1}\\right) p\\left(\\mathbf{x}_{1} | \\mathbf{z}_{1}\\right)=\\prod_{k=1}&#94;{K}\\left\\{\\pi_{k} p\\left(\\mathbf{x}_{1} | \\boldsymbol{\\phi}_{k}\\right)\\right\\}&#94;{z_{1 k}} & (13.37)\\\\ \\beta\\left(\\mathbf{z}_{n}\\right) &=\\sum_{\\mathbf{z}_{n+1}} \\beta\\left(\\mathbf{z}_{n+1}\\right) p\\left(\\mathbf{x}_{n+1} | \\mathbf{z}_{n+1}\\right) p\\left(\\mathbf{z}_{n+1} | \\mathbf{z}_{n}\\right) & (13.38)\\\\ \\beta\\left(\\mathbf{z}_{N}\\right) &= \\frac{\\gamma\\left(\\mathbf{z}_N\\right)p\\left(\\mathbf{X}\\right)}{\\alpha\\left(\\mathbf{z}_N\\right)} = \\frac{p\\left(\\mathbf{z}_{N} | \\mathbf{X}\\right)p(\\mathbf{X})}{p\\left(\\mathbf{X}, \\mathbf{z}_{N}\\right)} = 1 & (13.30) \\end{aligned} $$ \\(\\widehat{\\alpha},\\widehat{\\beta}\\) $$ \\begin{aligned} c_{n} \\widehat{\\alpha}\\left(\\mathbf{z}_{n}\\right) &=p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) \\sum_{\\mathbf{z}_{n-1}} \\widehat{\\alpha}\\left(\\mathbf{z}_{n-1}\\right) p\\left(\\mathbf{z}_{n} | \\mathbf{z}_{n-1}\\right) & (13.58)\\\\ \\widehat{\\alpha}\\left(\\mathbf{z}_{1}\\right) &=p\\left(\\mathbf{z}_{1}| \\mathbf{x}_{1}\\right)=\\frac{p\\left(\\mathbf{z}_{1}\\right)p\\left(\\mathbf{x}_{1}| \\mathbf{z}_{1}\\right)}{p\\left(\\mathbf{x}_{1}\\right)} = \\frac{\\prod_{k=1}&#94;{K}\\left\\{\\pi_{k} p\\left(\\mathbf{x}_{1} | \\boldsymbol{\\phi}_{k}\\right)\\right\\}&#94;{z_{1 k}}}{p\\left(\\mathbf{x}_{1}\\right)}\\\\ c_{n+1} \\widehat{\\beta}\\left(\\mathbf{z}_{n}\\right) &=\\sum_{\\mathbf{z}_{n+1}} \\widehat{\\beta}\\left(\\mathbf{z}_{n+1}\\right) p\\left(\\mathbf{x}_{n+1} | \\mathbf{z}_{n+1}\\right) p\\left(\\mathbf{z}_{n+1} | \\mathbf{z}_{n}\\right) &(13.62)\\\\ \\widehat{\\beta}\\left(\\mathbf{z}_{N}\\right) &= \\frac{\\gamma\\left(\\mathbf{z}_N\\right)}{\\widehat{\\alpha}\\left(\\mathbf{z}_N\\right)} = \\frac{p\\left(\\mathbf{z}_{N} | \\mathbf{X}\\right)}{p\\left(\\mathbf{z}_{n} | \\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n}\\right)} = 1 & (13.30) \\end{aligned} $$ なお、ここで \\((13.58)\\) でどのようにして \\(c_n\\) を求めるかですが、 $$ \\begin{aligned} \\mathrm{R.H.S}\\ (13.58) &= p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) \\sum_{\\mathbf{z}_{n-1}} \\widehat{\\alpha}\\left(\\mathbf{z}_{n-1}\\right) p\\left(\\mathbf{z}_{n} | \\mathbf{z}_{n-1}\\right)\\\\ &= p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) \\sum_{\\mathbf{z}_{n-1}}p\\left(\\mathbf{z}_{n-1}|\\mathbf{x}_1,\\ldots,\\mathbf{x}_{n-1}\\right) p\\left(\\mathbf{z}_{n} | \\mathbf{z}_{n-1}\\right)\\\\ &= p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) p\\left(\\mathbf{z}_{n}|\\mathbf{x}_1,\\ldots,\\mathbf{x}_{n-1}\\right)\\\\ &=p\\left(\\mathbf{x}_{n},\\mathbf{z}_{n}|\\mathbf{x}_1,\\ldots,\\mathbf{x}_{n-1}\\right) \\end{aligned} $$ となることから、 \\((13.58)\\) の右辺を \\(\\mathbf{z}_n\\) について周辺化すれば、 $$\\sum_{\\mathbf{z}_{n}}p\\left(\\mathbf{x}_{n},\\mathbf{z}_{n}|\\mathbf{x}_1,\\ldots,\\mathbf{x}_{n-1}\\right) = p\\left(\\mathbf{x}_{n}|\\mathbf{x}_1,\\ldots,\\mathbf{x}_{n-1}\\right) = c_n$$ となるので、 \\(c_n\\) が求められることがわかります。 おまけ（尤度関数） 尤度関数はスケーリング係数 \\(c\\) を用いるだけで簡単に求められることがわかります。 $$p(\\mathbf{X})=\\prod_{n=1}&#94;{N} c_{n}\\qquad (13.63)$$ if (!document.getElementById('mathjaxscript_pelican_#%@#$@#')) { var align = \"center\", indent = \"0em\", linebreak = \"false\"; if (false) { align = (screen.width < 768) ? \"left\" : align; indent = (screen.width < 768) ? \"0em\" : indent; linebreak = (screen.width < 768) ? 'true' : linebreak; } var mathjaxscript = document.createElement('script'); mathjaxscript.id = 'mathjaxscript_pelican_#%@#$@#'; mathjaxscript.type = 'text/javascript'; mathjaxscript.src = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.3/latest.js?config=TeX-AMS-MML_HTMLorMML'; var configscript = document.createElement('script'); configscript.type = 'text/x-mathjax-config'; configscript[(window.opera ? \"innerHTML\" : \"text\")] = \"MathJax.Hub.Config({\" + \" config: ['MMLorHTML.js'],\" + \" TeX: { extensions: ['AMSmath.js','AMSsymbols.js','noErrors.js','noUndefined.js'], equationNumbers: { autoNumber: 'none' } },\" + \" jax: ['input/TeX','input/MathML','output/HTML-CSS'],\" + \" extensions: ['tex2jax.js','mml2jax.js','MathMenu.js','MathZoom.js'],\" + \" displayAlign: '\"+ align +\"',\" + \" displayIndent: '\"+ indent +\"',\" + \" showMathMenu: true,\" + \" messageStyle: 'normal',\" + \" tex2jax: { \" + \" inlineMath: [ ['\\\\\\\\(','\\\\\\\\)'] ], \" + \" displayMath: [ ['$$','$$'] ],\" + \" processEscapes: true,\" + \" preview: 'TeX',\" + \" }, \" + \" 'HTML-CSS': { \" + \" fonts: [['STIX', 'TeX']],\" + \" styles: { '.MathJax_Display, .MathJax .mo, .MathJax .mi, .MathJax .mn': {color: 'inherit ! important'} },\" + \" linebreaks: { automatic: \"+ linebreak +\", width: '90% container' },\" + \" }, \" + \"}); \" + \"if ('default' !== 'default') {\" + \"MathJax.Hub.Register.StartupHook('HTML-CSS Jax Ready',function () {\" + \"var VARIANT = MathJax.OutputJax['HTML-CSS'].FONTDATA.VARIANT;\" + \"VARIANT['normal'].fonts.unshift('MathJax_default');\" + \"VARIANT['bold'].fonts.unshift('MathJax_default-bold');\" + \"VARIANT['italic'].fonts.unshift('MathJax_default-italic');\" + \"VARIANT['-tex-mathit'].fonts.unshift('MathJax_default-italic');\" + \"});\" + \"MathJax.Hub.Register.StartupHook('SVG Jax Ready',function () {\" + \"var VARIANT = MathJax.OutputJax.SVG.FONTDATA.VARIANT;\" + \"VARIANT['normal'].fonts.unshift('MathJax_default');\" + \"VARIANT['bold'].fonts.unshift('MathJax_default-bold');\" + \"VARIANT['italic'].fonts.unshift('MathJax_default-italic');\" + \"VARIANT['-tex-mathit'].fonts.unshift('MathJax_default-italic');\" + \"});\" + \"}\"; (document.body || document.getElementsByTagName('head')[0]).appendChild(configscript); (document.body || document.getElementsByTagName('head')[0]).appendChild(mathjaxscript); }","tags":"情報基礎実験","url":"https://iwasakishuto.github.io/University/3A/情報基礎実験-3.html","loc":"https://iwasakishuto.github.io/University/3A/情報基礎実験-3.html"},{"title":"HMMの最適化の計算過程","text":"ここでは、 HMMの最適化 で省略した計算過程について記述します。 \\(Q\\left(\\boldsymbol{\\theta}, \\boldsymbol{\\theta}&#94;{\\text {old }}\\right)\\) Maximization step $$ \\begin{aligned} Q\\left(\\boldsymbol{\\theta}, \\boldsymbol{\\theta}&#94;{\\mathrm{old}}\\right)=& \\sum_{k=1}&#94;{K} \\gamma\\left(z_{1 k}\\right) \\ln \\pi_{k}+\\sum_{n=2}&#94;{N} \\sum_{j=1}&#94;{K} \\sum_{k=1}&#94;{K} \\xi\\left(z_{n-1, j}, z_{n k}\\right) \\ln A_{j k} \\\\ &+\\sum_{n=1}&#94;{N} \\sum_{k=1}&#94;{K}\\gamma\\left(z_{n k}\\right) \\sum_{i=1}&#94;Dx_{ni} \\ln \\phi_{i k} \\end{aligned}\\qquad (13.17) $$ \\(\\boldsymbol{\\pi}\\) $$ \\begin{aligned} L\\left(\\boldsymbol{\\pi},\\boldsymbol{\\theta}, \\boldsymbol{\\theta}&#94;{\\text { old }}\\right) &= \\sum_{k=1}&#94;{K} \\gamma\\left(z_{1 k}\\right) \\ln \\pi_{k} + \\lambda_1\\left(\\sum_{k=1}&#94;K\\pi_k - 1\\right)\\\\ \\frac{\\partial L}{\\partial \\pi_k}&=\\frac{\\gamma\\left(z_{1 k}\\right)}{\\pi_k} + \\lambda_1 = 0\\quad \\therefore \\pi_k = -\\frac{\\gamma\\left(z_{1 k}\\right)}{\\lambda_1}\\\\ \\frac{\\partial L}{\\partial \\lambda_1}&=\\sum_{k=1}&#94;K\\pi_k - 1 = 0\\quad \\therefore\\lambda_1 = -\\sum_{k=1}&#94;K\\gamma\\left(z_{1 k}\\right)\\\\ \\therefore\\pi_k&#94;{\\star} &= \\frac{\\gamma\\left(z_{1 k}\\right)}{\\sum_{j=1}&#94;{K} \\gamma\\left(z_{1 j}\\right)} \\end{aligned} $$ \\(\\boldsymbol{A}\\) $$\\begin{aligned} L\\left(\\boldsymbol{A},\\boldsymbol{\\theta}, \\boldsymbol{\\theta}&#94;{\\text { old }}\\right) & = \\sum_{n=2}&#94;{N} \\sum_{j=1}&#94;{K} \\sum_{k=1}&#94;{K} \\xi\\left(z_{n-1, j}, z_{n k}\\right) \\ln A_{j k} + \\sum_{j=1}&#94;K\\lambda_{2,j}\\left(\\sum_{k=1}&#94;K A_{jk} - 1\\right)\\\\ \\frac{\\partial L}{\\partial A_{jk}} &= \\sum_{n=2}&#94;{N}\\frac{\\xi\\left(z_{n-1, j}, z_{n k}\\right)}{A_{jk}} + \\lambda_{2,j} = 0\\quad \\therefore A_{jk} = -\\frac{\\sum_{n=2}&#94;N\\xi\\left(z_{n-1, j}, z_{n k}\\right)}{\\lambda_{2,j}}\\\\ \\frac{\\partial L}{\\partial \\lambda_{2,j}} &= \\sum_{k=1}&#94;K A_{jk} - 1 = 0\\quad \\therefore \\lambda_{2,j} = -\\sum_{k=1}&#94;K\\sum_{n=2}&#94;N\\xi\\left(z_{n-1, j}, z_{n k}\\right)\\\\ \\therefore A_{jk}&#94;{\\star}&=\\frac{\\sum_{n=2}&#94;{N} \\xi\\left(z_{n-1, j}, z_{n k}\\right)}{\\sum_{l=1}&#94;{K} \\sum_{n=2}&#94;{N} \\xi\\left(z_{n-1, j}, z_{n l}\\right)} \\end{aligned}$$ \\(\\boldsymbol{\\phi}\\) $$\\begin{aligned} L\\left(\\boldsymbol{\\phi},\\boldsymbol{\\theta}, \\boldsymbol{\\theta}&#94;{\\text { old }}\\right) & = \\sum_{n=1}&#94;{N} \\sum_{k=1}&#94;{K}\\gamma\\left(z_{n k}\\right) \\sum_{i=1}&#94;Dx_{ni} \\ln \\phi_{i k} + \\sum_{k=1}&#94;K\\lambda_{3,k}\\left( \\sum_{i=1}&#94;D\\phi_{ik} - 1\\right)\\\\ \\frac{\\partial L}{\\partial\\phi_{ik}} &= \\frac{\\sum_{n=1}&#94;N\\gamma\\left(z_{n k}\\right)x_{ni}}{\\phi_{ik}} + \\lambda_{3,k}=0\\quad \\therefore \\phi_{ik} = -\\frac{\\sum_{n=1}&#94;N\\gamma\\left(z_{n k}\\right)x_{ni}}{\\lambda_{3,k}}\\\\ \\frac{\\partial L}{\\partial\\lambda_{3,k}} &= \\sum_{i=1}&#94;D\\phi_{ik} - 1 = 0 \\quad \\therefore \\lambda_{3,k} = -\\sum_{i=1}&#94;D\\sum_{n=1}&#94;N\\gamma\\left(z_{n k}\\right)x_{ni} = -\\sum_{i=1}&#94;D\\gamma\\left(z_{n k}\\right)\\\\ \\therefore \\phi_{ik}&#94;{\\star} &= \\frac{\\sum_{n=1}&#94;N\\gamma\\left(z_{n k}\\right)x_{ni}}{\\sum_{i=1}&#94;D\\gamma\\left(z_{n k}\\right)} \\end{aligned}$$ Expectation step $$ \\begin{aligned} \\gamma\\left(\\mathbf{z}_{n}\\right) &=p\\left(\\mathbf{z}_{n} | \\mathbf{X}, \\boldsymbol{\\theta}&#94;{\\text {old }}\\right) &(13.13)\\\\ \\xi\\left(\\mathbf{z}_{n-1}, \\mathbf{z}_{n}\\right) &=p\\left(\\mathbf{z}_{n-1}, \\mathbf{z}_{n} | \\mathbf{X}, \\boldsymbol{\\theta}&#94;{\\text {old }}\\right) &(13.14)\\\\ \\alpha\\left(\\mathbf{z}_{n}\\right) & \\equiv p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n}, \\mathbf{z}_{n}\\right) & (13.34)\\\\ \\beta\\left(\\mathbf{z}_{n}\\right) & \\equiv p\\left(\\mathbf{x}_{n+1}, \\ldots, \\mathbf{x}_{N} | \\mathbf{z}_{n}\\right) & (13.35)\\\\ \\end{aligned} $$ forward-algorithm \\((\\alpha)\\) $$\\begin{aligned} \\alpha\\left(\\mathbf{z}_{n}\\right) & =p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n}, \\mathbf{z}_{n}\\right) \\\\ & =p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) p\\left(\\mathbf{z}_{n}\\right) \\\\ & =p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n-1} | \\mathbf{z}_{n}\\right) p\\left(\\mathbf{z}_{n}\\right) \\quad (\\because \\text{conditional independence})\\\\ & =p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n-1}, \\mathbf{z}_{n}\\right) \\\\ & =p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) \\sum_{\\mathbf{z}_{n-1}} p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n-1}, \\mathbf{z}_{n-1}, \\mathbf{z}_{n}\\right) \\quad (\\because \\text{demarginalization})\\\\ & =p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) \\sum_{\\mathbf{z}_{n-1}} p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n-1}, \\mathbf{z}_{n} | \\mathbf{z}_{n-1}\\right) p\\left(\\mathbf{z}_{n-1}\\right) \\\\ & =p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) \\sum_{\\mathbf{z}_{n-1}} p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n-1} | \\mathbf{z}_{n-1}\\right) p\\left(\\mathbf{z}_{n} | \\mathbf{z}_{n-1}\\right) p\\left(\\mathbf{z}_{n-1}\\right) \\quad (\\because \\text{conditional independence})\\\\ & =p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) \\sum_{\\mathbf{z}_{n-1}} p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n-1}, \\mathbf{z}_{n-1}\\right) p\\left(\\mathbf{z}_{n} | \\mathbf{z}_{n-1}\\right) \\\\ & =p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) \\sum_{\\mathbf{z}_{n-1}} \\alpha(\\mathbf{z}_{n-1})p(\\mathbf{z} | \\mathbf{z}_{n-1})\\\\ \\end{aligned}$$ backward-algorithm \\((\\beta)\\) $$\\begin{aligned} \\beta\\left(\\mathbf{z}_{n}\\right) & =p\\left(\\mathbf{x}_{n+1}, \\ldots, \\mathbf{x}_{N} | \\mathbf{z}_{n}\\right) \\\\ & =\\sum_{\\mathbf{z}_{n+1}} p\\left(\\mathbf{x}_{n+1}, \\ldots, \\mathbf{x}_{N}, \\mathbf{z}_{n+1} | \\mathbf{z}_{n}\\right) \\quad (\\because \\text{demarginalization})\\\\ & =\\sum_{\\mathbf{z}_{n+1}} p\\left(\\mathbf{x}_{n+1}, \\ldots, \\mathbf{x}_{N} | \\mathbf{z}_{n}, \\mathbf{z}_{n+1}\\right) p\\left(\\mathbf{z}_{n+1} | \\mathbf{z}_{n}\\right) \\\\ & =\\sum_{\\mathbf{z}_{n+1}} p\\left(\\mathbf{x}_{n+1}, \\ldots, \\mathbf{x}_{N} | \\mathbf{z}_{n+1}\\right) p\\left(\\mathbf{z}_{n+1} | \\mathbf{z}_{n}\\right) \\quad (\\because \\text{conditional independence})\\\\ & =\\sum_{\\mathbf{z}_{n+1}} p\\left(\\mathbf{x}_{n+2}, \\ldots, \\mathbf{x}_{N} | \\mathbf{z}_{n+1}\\right) p\\left(\\mathbf{x}_{n+1} | \\mathbf{z}_{n+1}\\right) p\\left(\\mathbf{z}_{n+1} | \\mathbf{z}_{n}\\right) \\quad (\\because \\text{conditional independence})\\\\ & =\\sum_{\\mathbf{z}_{n+1}} \\beta(\\mathbf{z}_{n+1})p(\\mathbf{x}_{n+1}|\\mathbf{z}_{n+1})p(\\mathbf{z}_{n+1}|\\mathbf{z}_n)\\qquad (13.38) \\end{aligned}$$ \\(\\gamma,\\xi\\) $$ \\begin{aligned} \\gamma\\left(\\mathbf{z}_{n}\\right) &=p\\left(\\mathbf{z}_{n} | \\mathbf{X}\\right)=\\frac{p\\left(\\mathbf{X} | \\mathbf{z}_{n}\\right) p\\left(\\mathbf{z}_{n}\\right)}{p(\\mathbf{X})}\\quad (\\because \\text{Bayes' theorem}) &(13.32)\\\\ &=\\frac{p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n}, \\mathbf{z}_{n}\\right) p\\left(\\mathbf{x}_{n+1}, \\ldots, \\mathbf{x}_{N} | \\mathbf{z}_{n}\\right)}{p(\\mathbf{X})} \\quad (\\because \\text{conditional independence})\\\\ &= \\frac{\\alpha\\left(\\mathbf{z}_{n}\\right) \\beta\\left(\\mathbf{z}_{n}\\right)}{p(\\mathbf{X})} & (13.33) \\end{aligned} $$ $$\\begin{aligned} \\xi\\left(\\mathbf{z}_{n-1}, \\mathbf{z}_{n}\\right) &=p\\left(\\mathbf{z}_{n-1}, \\mathbf{z}_{n} | \\mathbf{X}\\right) \\\\ &=\\frac{p(\\mathbf{X} | \\mathbf{z}_{n-1}, \\mathbf{z}_{n}) p\\left(\\mathbf{z}_{n-1}, \\mathbf{z}_{n}\\right)}{p(\\mathbf{X})} \\quad (\\because \\text{Bayes' theorem})\\\\ &=\\frac{p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n-1} | \\mathbf{z}_{n-1}\\right) p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) p\\left(\\mathbf{x}_{n+1}, \\ldots, \\mathbf{x}_{N} | \\mathbf{z}_{n}\\right) p\\left(\\mathbf{z}_{n} | \\mathbf{z}_{n-1}\\right) p\\left(\\mathbf{z}_{n-1}\\right)}{p(\\mathbf{X})} \\quad (\\because \\text{conditional independence})\\\\ &=\\frac{\\alpha\\left(\\mathbf{z}_{n-1}\\right) p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) p\\left(\\mathbf{z}_{n} | \\mathbf{z}_{n-1}\\right) \\beta\\left(\\mathbf{x}_{n}\\right)}{p(\\mathbf{X})}\\qquad (13.43) \\end{aligned}$$ if (!document.getElementById('mathjaxscript_pelican_#%@#$@#')) { var align = \"center\", indent = \"0em\", linebreak = \"false\"; if (false) { align = (screen.width < 768) ? \"left\" : align; indent = (screen.width < 768) ? \"0em\" : indent; linebreak = (screen.width < 768) ? 'true' : linebreak; } var mathjaxscript = document.createElement('script'); mathjaxscript.id = 'mathjaxscript_pelican_#%@#$@#'; mathjaxscript.type = 'text/javascript'; mathjaxscript.src = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.3/latest.js?config=TeX-AMS-MML_HTMLorMML'; var configscript = document.createElement('script'); configscript.type = 'text/x-mathjax-config'; configscript[(window.opera ? \"innerHTML\" : \"text\")] = \"MathJax.Hub.Config({\" + \" config: ['MMLorHTML.js'],\" + \" TeX: { extensions: ['AMSmath.js','AMSsymbols.js','noErrors.js','noUndefined.js'], equationNumbers: { autoNumber: 'none' } },\" + \" jax: ['input/TeX','input/MathML','output/HTML-CSS'],\" + \" extensions: ['tex2jax.js','mml2jax.js','MathMenu.js','MathZoom.js'],\" + \" displayAlign: '\"+ align +\"',\" + \" displayIndent: '\"+ indent +\"',\" + \" showMathMenu: true,\" + \" messageStyle: 'normal',\" + \" tex2jax: { \" + \" inlineMath: [ ['\\\\\\\\(','\\\\\\\\)'] ], \" + \" displayMath: [ ['$$','$$'] ],\" + \" processEscapes: true,\" + \" preview: 'TeX',\" + \" }, \" + \" 'HTML-CSS': { \" + \" fonts: [['STIX', 'TeX']],\" + \" styles: { '.MathJax_Display, .MathJax .mo, .MathJax .mi, .MathJax .mn': {color: 'inherit ! important'} },\" + \" linebreaks: { automatic: \"+ linebreak +\", width: '90% container' },\" + \" }, \" + \"}); \" + \"if ('default' !== 'default') {\" + \"MathJax.Hub.Register.StartupHook('HTML-CSS Jax Ready',function () {\" + \"var VARIANT = MathJax.OutputJax['HTML-CSS'].FONTDATA.VARIANT;\" + \"VARIANT['normal'].fonts.unshift('MathJax_default');\" + \"VARIANT['bold'].fonts.unshift('MathJax_default-bold');\" + \"VARIANT['italic'].fonts.unshift('MathJax_default-italic');\" + \"VARIANT['-tex-mathit'].fonts.unshift('MathJax_default-italic');\" + \"});\" + \"MathJax.Hub.Register.StartupHook('SVG Jax Ready',function () {\" + \"var VARIANT = MathJax.OutputJax.SVG.FONTDATA.VARIANT;\" + \"VARIANT['normal'].fonts.unshift('MathJax_default');\" + \"VARIANT['bold'].fonts.unshift('MathJax_default-bold');\" + \"VARIANT['italic'].fonts.unshift('MathJax_default-italic');\" + \"VARIANT['-tex-mathit'].fonts.unshift('MathJax_default-italic');\" + \"});\" + \"}\"; (document.body || document.getElementsByTagName('head')[0]).appendChild(configscript); (document.body || document.getElementsByTagName('head')[0]).appendChild(mathjaxscript); }","tags":"情報基礎実験","url":"https://iwasakishuto.github.io/University/3A/情報基礎実験-2.html","loc":"https://iwasakishuto.github.io/University/3A/情報基礎実験-2.html"},{"title":"HMMの最適化","text":"隠れマルコフモデルは、一般に以下の式で表されます。 $$ p(\\mathbf{X}, \\mathbf{Z} | \\boldsymbol{\\theta})=p\\left(\\mathbf{z}_{1} | \\boldsymbol{\\pi}\\right)\\left[\\prod_{n=2}&#94;{N} p\\left(\\mathbf{z}_{n} | \\mathbf{z}_{n-1}, \\mathbf{A}\\right)\\right] \\prod_{m=1}&#94;{N} p\\left(\\mathbf{x}_{m} | \\mathbf{z}_{m}, \\boldsymbol{\\phi}\\right)\\qquad (13.10) $$ \\(\\mathbf{X} = \\{\\mathbf{x}_1,\\ldots,\\mathbf{x}_N\\}\\) : \\(\\mathbf{Z} = \\{\\mathbf{z}_1,\\ldots,\\mathbf{z}_N\\}\\) \\(\\boldsymbol{\\theta}=\\{\\boldsymbol{\\pi}, \\mathbf{A}, \\boldsymbol{\\phi}\\}\\) initial state \\(\\pi_{k} \\equiv p\\left(z_{1 k}=1\\right)\\) \\(p\\left(\\mathbf{z}_{1} \\mid \\pi\\right)=\\prod_{k=1}&#94;{K} \\pi_{k}&#94;{z_{1 k}}\\) transition probability \\(A_{j k} \\equiv p\\left(z_{n k}=1\\mid z_{n-1, j}=1\\right)\\) \\(p\\left(\\mathbf{z}_{n} \\mid \\mathbf{z}_{n-1}, \\mathbf{A}\\right)=\\prod_{k=1}&#94;{K} \\prod_{j=1}&#94;{K} A_{j k}&#94;{z_{n-1, j} z_{n k}}\\) emission probability \\(\\phi_{i k}\\equiv p\\left(x_{n i}=1 \\mid z_{n k}=1\\right)\\) \\(p(\\mathbf{x}_n \\mid \\mathbf{z}_n, \\boldsymbol{\\phi})=\\prod_{i=1}&#94;{D} \\prod_{k=1}&#94;{K} \\phi_{i k}&#94;{x_{ni} z_{nk}}\\) ※ 一般に放出確率はどのような分布を考えることも可能ですが、今回は上記の離散多項分布を考えます。 尤度関数 ここで、データ集合 \\(\\mathbf{X}\\) が観測された際に、上記の同時分布を潜在変数 \\(\\mathbf{Z}\\) について周辺化することで、尤度関数は以下のように記述されます。 $$ p(\\mathbf{X} | \\boldsymbol{\\theta})=\\sum_{\\mathbf{Z}} p(\\mathbf{X}, \\mathbf{Z} | \\boldsymbol{\\theta})\\qquad (13.11) $$ しかし、この尤度関数は \\(n\\) について分解できない（ \\(\\mathbf{z}_n\\) ）ごとに和を取れないので、 条件付き独立 の性質を活かして 尤度関数の対数の期待値 を最大化する Baum-Welch algorithm (EM algorithm) を用います。 Baum-Welch (EM) パラメータ \\(\\boldsymbol{\\theta}&#94;{\\text {old }}\\) を用いて \\(p\\left(\\mathbf{Z} | \\mathbf{X}, \\boldsymbol{\\theta}&#94;{\\text {old }}\\right)\\) を最大化する。 対数尤度関数の期待値 \\(Q\\left(\\boldsymbol{\\theta}, \\boldsymbol{\\theta}&#94;{\\text {old }}\\right)\\) を求める。 \\(Q\\left(\\boldsymbol{\\theta}, \\boldsymbol{\\theta}&#94;{\\text {old }}\\right)\\) を最大化するパラメータに更新する。 \\(\\boldsymbol{\\theta}\\rightarrow\\boldsymbol{\\theta}&#94;{\\text {old }}\\) 1に戻る。 $$ Q\\left(\\boldsymbol{\\theta}, \\boldsymbol{\\theta}&#94;{\\text {old }}\\right)=\\sum_{\\mathbf{Z}} p\\left(\\mathbf{Z} | \\mathbf{X}, \\boldsymbol{\\theta}&#94;{\\text {old }}\\right) \\ln p(\\mathbf{X}, \\mathbf{Z} | \\boldsymbol{\\theta})\\qquad (13.12) $$ ここで、表記を簡単にするために、 γ 、 ξ を導入します。 $$ \\begin{aligned} \\gamma\\left(\\mathbf{z}_{n}\\right) &=p\\left(\\mathbf{z}_{n} | \\mathbf{X}, \\boldsymbol{\\theta}&#94;{\\text {old }}\\right) &(13.13)\\\\ \\xi\\left(\\mathbf{z}_{n-1}, \\mathbf{z}_{n}\\right) &=p\\left(\\mathbf{z}_{n-1}, \\mathbf{z}_{n} | \\mathbf{X}, \\boldsymbol{\\theta}&#94;{\\text {old }}\\right) &(13.14) \\end{aligned} $$ なお、潜在変数が離散なので、以下のように記述し直せます。（ \\(\\pi_{k},A_{j k},\\phi_{j k}\\) と同様。） $$ \\begin{aligned} \\gamma\\left(z_{n k}\\right) &=\\sum_{\\mathbf{z}} \\gamma(\\mathbf{z}) z_{n k} &(13.15)\\\\ \\xi\\left(z_{n-1, j}, z_{n k}\\right) &=\\sum_{\\mathbf{z}} \\gamma(\\mathbf{z}) z_{n-1, j} z_{n k} &(13.16)\\end{aligned} $$ これらを用いると、 \\(Q\\left(\\boldsymbol{\\theta}, \\boldsymbol{\\theta}&#94;{\\mathrm{old}}\\right)\\) が以下のように書き下せます。 （計算過程） $$ \\begin{aligned} Q\\left(\\boldsymbol{\\theta}, \\boldsymbol{\\theta}&#94;{\\mathrm{old}}\\right)=& \\sum_{k=1}&#94;{K} \\gamma\\left(z_{1 k}\\right) \\ln \\pi_{k}+\\sum_{n=2}&#94;{N} \\sum_{j=1}&#94;{K} \\sum_{k=1}&#94;{K} \\xi\\left(z_{n-1, j}, z_{n k}\\right) \\ln A_{j k} \\\\ &+\\sum_{n=1}&#94;{N} \\sum_{k=1}&#94;{K}\\gamma\\left(z_{n k}\\right) \\sum_{i=1}&#94;Dx_{ni} \\ln \\phi_{i k} \\end{aligned}\\qquad (13.17) $$ Maximization step ※ 実際の計算の順番からは前後しますが、先にM stepを説明します。 上記の \\(Q\\left(\\boldsymbol{\\theta}, \\boldsymbol{\\theta}&#94;{\\mathrm{old}}\\right)\\) を各パラメータ \\(\\boldsymbol{\\theta}\\) に関して最大化するのは（別ブロックに別れているから）簡単で、それぞれ 適当なラグランジュ乗数を導入する ことで、以下のように求まります。 （計算過程） $$ \\begin{aligned} \\pi_{k}&= \\frac{\\gamma\\left(z_{1 k}\\right)}{\\sum_{j=1}&#94;{K} \\gamma\\left(z_{1 j}\\right)} & (13.18)\\\\ A_{j k}&= \\frac{\\sum_{n=2}&#94;{N} \\xi\\left(z_{n-1, j}, z_{n k}\\right)}{\\sum_{l=1}&#94;{K} \\sum_{n=2}&#94;{N} \\xi\\left(z_{n-1, j}, z_{n l}\\right)} & (13.19)\\\\ \\phi_{i k}&=\\frac{\\sum_{n=1}&#94;{N} \\gamma\\left(z_{n k}\\right) x_{n i}}{\\sum_{n=1}&#94;{N} \\gamma\\left(z_{n k}\\right)} & (13.23) \\end{aligned} $$ Expectation step M step で必要となる \\(\\gamma,\\xi\\) は、 条件付き独立性 を用いることで効率的に計算することができます。 $$ \\begin{aligned} \\gamma\\left(\\mathbf{z}_{n}\\right) &=p\\left(\\mathbf{z}_{n} | \\mathbf{X}\\right)=\\frac{p\\left(\\mathbf{X} | \\mathbf{z}_{n}\\right) p\\left(\\mathbf{z}_{n}\\right)}{p(\\mathbf{X})}\\quad (\\because \\text{Bayes' theorem}) &(13.32)\\\\ &=\\frac{p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n}, \\mathbf{z}_{n}\\right) p\\left(\\mathbf{x}_{n+1}, \\ldots, \\mathbf{x}_{N} | \\mathbf{z}_{n}\\right)}{p(\\mathbf{X})} \\quad (\\because \\text{conditional independence})\\\\ &= \\frac{\\alpha\\left(\\mathbf{z}_{n}\\right) \\beta\\left(\\mathbf{z}_{n}\\right)}{p(\\mathbf{X})} & (13.33) \\end{aligned} $$ $$ \\begin{aligned} \\alpha\\left(\\mathbf{z}_{n}\\right) & \\equiv p\\left(\\mathbf{x}_{1}, \\ldots, \\mathbf{x}_{n}, \\mathbf{z}_{n}\\right) & (13.34)\\\\ \\beta\\left(\\mathbf{z}_{n}\\right) & \\equiv p\\left(\\mathbf{x}_{n+1}, \\ldots, \\mathbf{x}_{N} | \\mathbf{z}_{n}\\right) & (13.35) \\end{aligned} $$ 条件付き独立性を用いてそれぞれ変形すると、以下の再帰式を導くことができます。 （計算過程） $$ \\begin{aligned} \\alpha\\left(\\mathbf{z}_{n}\\right) &=p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) \\sum_{\\mathbf{z}_{n-1}} \\alpha\\left(\\mathbf{z}_{n-1}\\right) p\\left(\\mathbf{z}_{n} | \\mathbf{z}_{n-1}\\right) & (13.36)\\\\ \\alpha\\left(\\mathbf{z}_{1}\\right) &=p\\left(\\mathbf{x}_{1}, \\mathbf{z}_{1}\\right)=p\\left(\\mathbf{z}_{1}\\right) p\\left(\\mathbf{x}_{1} | \\mathbf{z}_{1}\\right)=\\prod_{k=1}&#94;{K}\\left\\{\\pi_{k} p\\left(\\mathbf{x}_{1} | \\boldsymbol{\\phi}_{k}\\right)\\right\\}&#94;{z_{1 k}} & (13.37)\\\\ \\beta\\left(\\mathbf{z}_{n}\\right) &=\\sum_{\\mathbf{z}_{n+1}} \\beta\\left(\\mathbf{z}_{n+1}\\right) p\\left(\\mathbf{x}_{n+1} | \\mathbf{z}_{n+1}\\right) p\\left(\\mathbf{z}_{n+1} | \\mathbf{z}_{n}\\right) & (13.38)\\\\ \\beta\\left(\\mathbf{z}_{N}\\right) &= \\frac{\\gamma\\left(\\mathbf{z}_N\\right)p\\left(\\mathbf{X}\\right)}{\\alpha\\left(\\mathbf{z}_N\\right)} = \\frac{p\\left(\\mathbf{z}_{N} | \\mathbf{X}\\right)p(\\mathbf{X})}{p\\left(\\mathbf{X}, \\mathbf{z}_{N}\\right)} = 1 & (13.30) \\end{aligned} $$ ※ なお、ここで \\(\\alpha\\) の再起式を forward-algorithm 、 \\(\\beta\\) の再起式を backward-algorithm と呼ぶことがあります。 また、これらを用いて \\(\\xi\\) を記述することもでき、以下のように表せます。 （計算過程） $$\\xi\\left(\\mathbf{z}_{n-1}, \\mathbf{z}_{n}\\right)=\\frac{\\alpha\\left(\\mathbf{z}_{n-1}\\right) p\\left(\\mathbf{x}_{n} | \\mathbf{z}_{n}\\right) p\\left(\\mathbf{z}_{n} | \\mathbf{z}_{n-1}\\right) \\beta\\left(\\mathbf{z}_{n}\\right)}{p(\\mathbf{X})}\\qquad (13.43)$$ 以上でBaum-Welchに必要な計算式が求まりました。 おまけ（尤度関数） 尤度関数は、アルゴリズムの停止条件に用いられるなど、値を求めることが非常に有用です。 求める際は、 $$ \\gamma\\left(\\mathbf{z}_{n}\\right)= \\frac{\\alpha\\left(\\mathbf{z}_{n}\\right) \\beta\\left(\\mathbf{z}_{n}\\right)}{p(\\mathbf{X})} \\qquad (13.33) $$ の両辺を \\(\\mathbf{z}_n\\) について周辺化すれば、左辺は $$\\sum_{\\mathbf{z}_{n}} \\gamma\\left(\\mathbf{z}_{n}\\right) = \\sum_{\\mathbf{z}_{n}}p\\left(\\mathbf{z}_{n} | \\mathbf{X}\\right) = 1$$ となることが明らかなので、以下のように求まります。 $$p(\\mathbf{X})=\\sum_{\\mathbf{z}_{n}} \\alpha\\left(\\mathbf{z}_{n}\\right) \\beta\\left(\\mathbf{z}_{n}\\right)\\qquad (13.41)$$ また、上記の式は任意の \\(n\\) について成立するので、 \\(n=N\\) の場合を考えれば \\(\\alpha\\) のみを用いて求めることができます。 $$p(\\mathbf{X})=\\sum_{\\mathbf{z}_{N}} \\alpha\\left(\\mathbf{z}_{N}\\right)\\qquad (13.42)$$ if (!document.getElementById('mathjaxscript_pelican_#%@#$@#')) { var align = \"center\", indent = \"0em\", linebreak = \"false\"; if (false) { align = (screen.width < 768) ? \"left\" : align; indent = (screen.width < 768) ? \"0em\" : indent; linebreak = (screen.width < 768) ? 'true' : linebreak; } var mathjaxscript = document.createElement('script'); mathjaxscript.id = 'mathjaxscript_pelican_#%@#$@#'; mathjaxscript.type = 'text/javascript'; mathjaxscript.src = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.3/latest.js?config=TeX-AMS-MML_HTMLorMML'; var configscript = document.createElement('script'); configscript.type = 'text/x-mathjax-config'; configscript[(window.opera ? \"innerHTML\" : \"text\")] = \"MathJax.Hub.Config({\" + \" config: ['MMLorHTML.js'],\" + \" TeX: { extensions: ['AMSmath.js','AMSsymbols.js','noErrors.js','noUndefined.js'], equationNumbers: { autoNumber: 'none' } },\" + \" jax: ['input/TeX','input/MathML','output/HTML-CSS'],\" + \" extensions: ['tex2jax.js','mml2jax.js','MathMenu.js','MathZoom.js'],\" + \" displayAlign: '\"+ align +\"',\" + \" displayIndent: '\"+ indent +\"',\" + \" showMathMenu: true,\" + \" messageStyle: 'normal',\" + \" tex2jax: { \" + \" inlineMath: [ ['\\\\\\\\(','\\\\\\\\)'] ], \" + \" displayMath: [ ['$$','$$'] ],\" + \" processEscapes: true,\" + \" preview: 'TeX',\" + \" }, \" + \" 'HTML-CSS': { \" + \" fonts: [['STIX', 'TeX']],\" + \" styles: { '.MathJax_Display, .MathJax .mo, .MathJax .mi, .MathJax .mn': {color: 'inherit ! important'} },\" + \" linebreaks: { automatic: \"+ linebreak +\", width: '90% container' },\" + \" }, \" + \"}); \" + \"if ('default' !== 'default') {\" + \"MathJax.Hub.Register.StartupHook('HTML-CSS Jax Ready',function () {\" + \"var VARIANT = MathJax.OutputJax['HTML-CSS'].FONTDATA.VARIANT;\" + \"VARIANT['normal'].fonts.unshift('MathJax_default');\" + \"VARIANT['bold'].fonts.unshift('MathJax_default-bold');\" + \"VARIANT['italic'].fonts.unshift('MathJax_default-italic');\" + \"VARIANT['-tex-mathit'].fonts.unshift('MathJax_default-italic');\" + \"});\" + \"MathJax.Hub.Register.StartupHook('SVG Jax Ready',function () {\" + \"var VARIANT = MathJax.OutputJax.SVG.FONTDATA.VARIANT;\" + \"VARIANT['normal'].fonts.unshift('MathJax_default');\" + \"VARIANT['bold'].fonts.unshift('MathJax_default-bold');\" + \"VARIANT['italic'].fonts.unshift('MathJax_default-italic');\" + \"VARIANT['-tex-mathit'].fonts.unshift('MathJax_default-italic');\" + \"});\" + \"}\"; (document.body || document.getElementsByTagName('head')[0]).appendChild(configscript); (document.body || document.getElementsByTagName('head')[0]).appendChild(mathjaxscript); }","tags":"情報基礎実験","url":"https://iwasakishuto.github.io/University/3A/情報基礎実験-1.html","loc":"https://iwasakishuto.github.io/University/3A/情報基礎実験-1.html"}]};