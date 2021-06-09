### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 77f95570-c78c-11eb-2d0c-d12e6e2ad622
begin
	using Pkg
	Pkg.activate(pwd());
	Pkg.instantiate();
	# Pkg.add(["PlutoUI","LaTeXStrings","Plots","Images"])
	
	using PlutoUI, LaTeXStrings, Plots
end

# ╔═╡ d4b63e91-431d-4097-8a11-0e3eb404b91b
PlutoUI.TableOfContents(depth = 10)

# ╔═╡ 5798e810-5238-4805-8f74-7448be1792b4
md"### resources"

# ╔═╡ e642d15f-7539-46a8-a302-3b04ccbfc6ef
md"### structure"

# ╔═╡ 47dc183e-78c0-49ce-a0a1-e0a417923bca
md"### function"

# ╔═╡ 587c84ce-bb7f-4f9b-9d5e-4f79ec4790fa
function returnperiod(T,P)
	return R = abs(T/log(1-P))
end

# ╔═╡ bb853a92-2571-4943-9f92-aa73f778ea1c
struct Earthquake
	code
	title
	label
	Ss # 震區短週期水平譜加速度係數
	S1 # 震區1秒週期水平譜加速度係數
	NA # 近斷層調整因子
	NV
	Fa # 短週期結構之工址放大係數
	Fv # 長週期結構之工址放大係數
	Ss_site # 近斷層區域工址短週期水平譜加速度係數
	S1_site # 近斷層區域工址1秒週期水平譜加速度係數
	T0
	T
	Sa
	T_Prob
	Prob
	R # 回歸期
	function Earthquake(code::String,Ss::N,S1::N,NA::N,NV::N,Fa::N,Fv::N) where {N<:Real}
		S_short = Ss*Fa*NA;
		S_1sec = S1*Fv*NV;
		T0 = S_1sec/S_short;
		T = collect(range(0.01,5T0,length = 1000));
		
		Sa = fill(NaN, 1000);
		Sa[T .<= 0.2T0] .= S_short*(0.4 .+ 3*T[T .<= 0.2T0]/T0);
		Sa[0.2T0 .< T .<= T0] .= S_short;
		Sa[1T0 .< T .<= 2.5T0] .= S_1sec ./ T[1T0 .< T .<= 2.5T0];
		Sa[T .> 2.5T0] .= 0.4*S_short;
		
		# about return period
		T_Prob = 50;
		if code == "M"
			label = "max";
			title = "最大考量地震";
			Prob = 0.02;
		elseif code == "D"
			label = "design";
			title = "設計地震";
			Prob = 0.1;
		else
			error(""" Unsupported code name. It should be either "M" or "D". """);
		end
		R = returnperiod(T_Prob,Prob);
		
		return new(code,title,label,Ss,S1,NA,NV,Fa,Fv,S_short,S_1sec,T0,T,Sa,T_Prob,Prob,R);
		
	end
end

# ╔═╡ ebb47979-43a5-4396-9b4b-d5d67fb76a2b
function plot(Eq::Earthquake; c...)
	p = Plots.plot(Eq.T, Eq.Sa; c...);
	return p
end

# ╔═╡ 52915591-d7c8-4cec-9b4f-3c0c54efeb5f
function plot!(p, Eq::Earthquake; c...)
	Plots.plot!(p, Eq.T, Eq.Sa; c...);
end

# ╔═╡ f438c670-5672-4a1e-87e0-7ca73abcd389
function show_returnperiod(Eq::Earthquake)
if Eq.code == "D"
	dscp1 = "容許塑性變形但韌性需求不得超過容許韌性容量";
elseif Eq.code == "M"
	dscp1 = "使用之韌性可達規定之韌性容量";
end
	
R_approx = Int64(round(Eq.R));
		
s = """
#### $(Eq.title)
- $(dscp1)
-  ``T = $(Eq.T_Prob)`` 年內超越機率``P> $(Eq.Prob)``； 即回歸期 ``R\\approx $R_approx`` 年 ( ``$(Eq.R)`` 年)
"""
Markdown.parse(s)
end

# ╔═╡ b127b37a-1bf9-4f7d-9dc4-451d162bc471
function show_parameter(Eq::Earthquake)
	Mk = Eq.code;
	str = """
	查表參數：
	- 震區短週期水平譜加速度係數: ``S_S^$(Mk) = $(Eq.Ss)``; 
	- 震區1秒週期水平譜加速度係數: ``S_1^$(Mk) = $(Eq.S1)``; 
	- 近斷層調整因子: ``N_A = $(Eq.NA)``; ``N_V = $(Eq.NV)``; 
	- 短週期結構之工址放大係數: ``Fa = $(Eq.Fa)``
	- 長週期結構之工址放大係數: ``Fv = $(Eq.Fv)`` 
	
	計算結果：
	- 近斷層區域工址短週期水平譜加速度係數: ``S_{$(Mk)S} = S_S^$(Mk) F_a N_A = $(Eq.Ss_site)`` 
	- 近斷層區域工址1秒週期水平譜加速度係數: ``S_{$(Mk)1} = S_1^$(Mk) F_v N_V = $(Eq.S1_site)`` 
	- ``T_0^$(Mk) = \\frac{S_{$(Mk)1}}{S_{$(Mk)S}} =  $(Eq.T0)``

	"""
	# Ss::N # 震區短週期水平譜加速度係數
	# S1::N # 震區1秒週期水平譜加速度係數
	# NA::N # 近斷層調整因子
	# NV::N
	# Fa::N # 短週期結構之工址放大係數
	# Fv::N # 長週期結構之工址放大係數
	# Ss_site # 近斷層區域工址短週期水平譜加速度係數
	# S1_site # 近斷層區域工址1秒週期水平譜加速度係數
	@eval @md_str $str
end

# ╔═╡ 752508a2-e1ed-42d3-85f6-a2ee9f765203
@bind refResh Button("Refresh")

# ╔═╡ 823f2e7a-5ade-43dc-b844-d8223f01bfa3
md"### Input parameter"

# ╔═╡ 60cc0a5a-1efa-46a3-aa87-e63a9b4f11bd
Eq_design = Earthquake("D",0.7,0.4,1.03,1.1,1.0,1.0)

# ╔═╡ afc2727e-8a44-40ad-a6a3-d2da3c98f8cb
Eq_max = Earthquake("M",0.9,0.5,1.03,1.15,1.0,1.0)

# ╔═╡ b81dcc3f-d566-492e-a0e0-1c339424156e
with_terminal() do
dump(Eq_design)
end

# ╔═╡ 76b0ce07-68f2-429f-9515-3188a8850ce2
with_terminal() do
dump(Eq_max)
end

# ╔═╡ 6b37da3c-8e16-45dc-b86b-37ea6acda46e
md"# 地震工程實務分析期末作業
107682001 吳宗羲"

# ╔═╡ b1911a8e-4d6e-4bd7-9aba-7219dc32bf77
md"""
### 基本資訊
#### 回歸期公式

$R = T/\log(1-P(Z>z))$
"""

# ╔═╡ d084e446-a1f7-41ea-a858-08ab24a85c7c
show_returnperiod(Eq_design)

# ╔═╡ a6a03966-be99-445a-ace3-6e3373ea9d00
show_returnperiod(Eq_max)

# ╔═╡ 2970a59e-2231-427b-aa6c-354aeb0083b1
md"### 參數"

# ╔═╡ 5169a46f-a454-4d37-86a3-2c91fbb3988a
md"""
#### 工址
![](https://github.com/okatsn/FinalHW_EarthquakeEngineering/raw/master/img/site_loc.png)

#### 震區水平譜加速度係數
  
![](https://github.com/okatsn/FinalHW_EarthquakeEngineering/raw/master/img/Fig_table_2-1.png)

#### 近斷層調整因子 (以車籠埔斷層為例)
![](https://github.com/okatsn/FinalHW_EarthquakeEngineering/raw/master/img/Fig_table_2-4-1.png)
"""

# ╔═╡ d82c1e80-60f0-4fb7-a709-aab2d948c8b9
md"#### 工址放大係數 (第一類地盤)"

# ╔═╡ 910cd27a-d11e-41bc-a202-372f1d8ebc0f
md"### 計算結果"

# ╔═╡ 852f3717-84ec-44f9-850d-b244f3b5a089
md"#### 設計地震"

# ╔═╡ c128d099-333a-461d-ab26-aa43a4bde5a0
begin
	refResh
	show_parameter(Eq_design)
end

# ╔═╡ 7fe0908d-9d38-4198-85cd-7849868bface
md"#### 最大參考地震"

# ╔═╡ 8c51288d-8d6e-49ba-918b-af9d9088a8c8
begin
	refResh
	show_parameter(Eq_max)
end

# ╔═╡ e4aa31b8-7a39-4a61-9ece-38c77755a6e7
let 
	p = plot(Eq_design,labels = "design");
	plot!(p, Eq_max, labels = "max", size=(600,300), xlabel = "T", ylabel = "Sa");
end

# ╔═╡ Cell order:
# ╠═d4b63e91-431d-4097-8a11-0e3eb404b91b
# ╠═77f95570-c78c-11eb-2d0c-d12e6e2ad622
# ╠═5798e810-5238-4805-8f74-7448be1792b4
# ╟─e642d15f-7539-46a8-a302-3b04ccbfc6ef
# ╠═bb853a92-2571-4943-9f92-aa73f778ea1c
# ╟─47dc183e-78c0-49ce-a0a1-e0a417923bca
# ╠═ebb47979-43a5-4396-9b4b-d5d67fb76a2b
# ╠═52915591-d7c8-4cec-9b4f-3c0c54efeb5f
# ╠═587c84ce-bb7f-4f9b-9d5e-4f79ec4790fa
# ╠═f438c670-5672-4a1e-87e0-7ca73abcd389
# ╠═b127b37a-1bf9-4f7d-9dc4-451d162bc471
# ╟─752508a2-e1ed-42d3-85f6-a2ee9f765203
# ╟─823f2e7a-5ade-43dc-b844-d8223f01bfa3
# ╠═60cc0a5a-1efa-46a3-aa87-e63a9b4f11bd
# ╠═afc2727e-8a44-40ad-a6a3-d2da3c98f8cb
# ╠═b81dcc3f-d566-492e-a0e0-1c339424156e
# ╠═76b0ce07-68f2-429f-9515-3188a8850ce2
# ╟─6b37da3c-8e16-45dc-b86b-37ea6acda46e
# ╟─b1911a8e-4d6e-4bd7-9aba-7219dc32bf77
# ╟─d084e446-a1f7-41ea-a858-08ab24a85c7c
# ╟─a6a03966-be99-445a-ace3-6e3373ea9d00
# ╟─2970a59e-2231-427b-aa6c-354aeb0083b1
# ╠═5169a46f-a454-4d37-86a3-2c91fbb3988a
# ╟─d82c1e80-60f0-4fb7-a709-aab2d948c8b9
# ╟─910cd27a-d11e-41bc-a202-372f1d8ebc0f
# ╟─852f3717-84ec-44f9-850d-b244f3b5a089
# ╟─c128d099-333a-461d-ab26-aa43a4bde5a0
# ╟─7fe0908d-9d38-4198-85cd-7849868bface
# ╟─8c51288d-8d6e-49ba-918b-af9d9088a8c8
# ╠═e4aa31b8-7a39-4a61-9ece-38c77755a6e7