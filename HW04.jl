### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ db61ecb1-6ed9-44d9-b8e4-b50622efbac1
using TestImages, ImageShow, Random, PlutoTest, Plots, FourierTools, Statistics, FFTW, IndexFunArrays, SpecialFunctions, PlutoUI, Noise, Optim, ForwardDiff, ImageIO, ImageShow

# ╔═╡ fbc3ec99-0c73-40e7-bd21-933bb5a52755
begin
	
	"""
	    simshow(arr; set_one=false, set_zero=false,
	                f=nothing, γ=1)
	Displays a real valued array . Brightness encodes magnitude.
	Works within Jupyter and Pluto.
	# Keyword args
	The transforms are applied in that order.
	* `set_zero=false` subtracts the minimum to set minimum to 1
	* `set_one=false` divides by the maximum to set maximum to 1
	* `f` applies an arbitrary function to the abs array
	* `γ` applies a gamma correction to the abs 
	* `cmap=:gray` applies a colormap provided by ColorSchemes.jl. If `cmap=:gray` simply `Colors.Gray` is used
	    and with different colormaps the result is an `Colors.RGB` element type
	"""
	function simshow(arr::AbstractArray{<:Real}; 
	                 set_one=true, set_zero=false,
	                 f = nothing,
	                 γ = 1,
	                 cmap=:gray)
	    arr = set_zero ? arr .- minimum(arr) : arr
	
	    if set_one
	        m = maximum(arr)
	        if !iszero(m)
	            arr = arr ./ maximum(arr)
	        end
	    end
	
	    arr = isnothing(f) ? arr : f(arr)
	
	    if !isone(γ)
	        arr = arr .^ γ
	    end
	
	
	    if cmap == :gray
	        Gray.(arr)
	    else
	        get(colorschemes[cmap], arr)
	    end
	end
	
	
	"""
	    simshow(arr)
	Displays a complex array. Color encodes phase, brightness encodes magnitude.
	Works within Jupyter and Pluto.
	# Keyword args
	The transforms are applied in that order.
	* `f` applies a function `f` to the array.
	* `absf` applies a function `absf` to the absolute of the array
	* `absγ` applies a gamma correction to the abs 
	"""
	function simshow(arr::AbstractArray{T};
	                 f=nothing,
	                 absγ=1,
	                 absf=nothing) where (T<:Complex)
	
	    if !isnothing(f)
	        arr = f(arr)
	    end
	
	    Tr = real(T)
	    # scale abs to 1
	    absarr = abs.(arr)
	    absarr ./= maximum(absarr)
	
	    if !isnothing(absf)
	        absarr .= absf(absarr)
	    end
	    
	    if !isone(absγ)
	        absarr .= absarr .^ absγ
	    end
	
	    angarr = angle.(arr) ./ Tr(2pi) * 360 
	
	    HSV.(angarr, one(Tr), absarr)
	end
	
	
	
	"""
	    simshow(arr::AbstractArray{<:Colors.ColorTypes.Colorant})
	If `simshow` receives an array which already contains color information, just display it.
	"""
	function simshow(arr::AbstractArray{<:Colors.ColorTypes.Colorant})
	    return arr
	end
end

# ╔═╡ 82726a23-060c-4c7d-91c1-e78a22fb98d6
md"# 1 Fourier Transforms

In the last homework we learnt that FFTW calculates the Fourier transform such that the DC frequency (zero frequency) is at the first entry of the array.

To shift the frequency to the center, we normally use `fftshift`.
However, another issue is, that the FFT also interpretes the center of the array at the first index position.
So for example, we would expect that the Fourier transform of the delta peak is a constant array! 

As you can see below, the Fourier transform is only a constant array if the delta peak is at the fist entry
"

# ╔═╡ ebe90686-c72f-4a42-a644-f06babe00f41
# FFT of delta peak results in constant array
fft([1.0, 0.0, 0.0])

# ╔═╡ 5053c100-5f27-4f2f-acbe-e12c9070feac
# no constant array since delta is not at the first entry
fft([0.0, 1.0, 0.0])

# ╔═╡ 361a7a34-d70b-4597-948c-80ed913a562b
# no constant array since delta is not at the first entry
fft([0.0, 0.0, 1.0])

# ╔═╡ 653a94c9-cb89-47df-8ba5-0ce683440b78
md"
The reason for this behaviour is that the second example is a shifted array of the first one. Hence, in Fourier space this corresponds to a phase ramp!

In the last homework we programmed `ffts(x) = fftshift(fft(x))` and `iffts(x) = ifft(ifftshift(x))`.

We now also introduce `my_ft` and `my_ift`.


The qualitative meaning of the different conventions for the Fourier transform of a signal with length `N` is:
* `fft`: center in real space is at `(1,)` and in Fourier space at `(1,)`
* `ffts`: center in real space is at `(1,)` and in Fourier space at `(N ÷ 2 + 1,)`
* `ft`: center in real space is at `(N ÷ 2 + 1,)` and in Fourier space at `(N ÷ 2 + 1,)`


"

# ╔═╡ b65e0745-82c6-48e1-a956-fa71583dad15
begin
	my_ft(x) = fftshift(fft(ifftshift(x)))
	my_ift(x) = fftshift(ifft(ifftshift(x)))
end

# ╔═╡ c53325bd-eeb6-4663-867d-44fd3ace273b
md"## Task 1.1
Try to change the following test below.
You are only allowed to insert `fftshift` and `ifftshift` statements.
Don't change the order of the Fourier transform calls.
"

# ╔═╡ b57917e8-3d2e-433b-84df-6129f760f954
begin
	arr_even = [1.0,2,3,4]
	arr_odd = [1.0,2,3,4,5]
end

# ╔═╡ 6e0e223d-9107-4121-880a-8ed5d7e5e9c2
md"The Test is broken but we fix it by changing the right hand side. Do that for the red tests below as well"

# ╔═╡ 5213f82b-5e06-4cee-a0d1-21f1c5cc2998
PlutoTest.@test ffts(arr_even) ≈ fft(arr_even)

# ╔═╡ 64103e19-ce79-48b7-bac9-b429f6d423c2
# this one is fixed
PlutoTest.@test ffts(arr_even) ≈ fftshift(fft(arr_even))

# ╔═╡ f380db27-8773-413b-978f-4496fc585ae3
md"### Task 1.1
Now try to fix always the right hand side accordingly!
"

# ╔═╡ 0464d972-fcaf-4f83-be7d-073d217e8f4c
# TODO
PlutoTest.@test ift(arr_odd) ≈ iffts(arr_odd)

# ╔═╡ 87fabbfe-e25b-467a-9fbd-ceec98ba4ed6
# TODO
PlutoTest.@test iffts(arr_odd) ≈ ifft(arr_odd)

# ╔═╡ c537e21b-5ed2-4cac-947a-e052474a2442
# TODO
PlutoTest.@test ffts(arr_odd) ≈ ft(arr_odd)

# ╔═╡ 3093ae8d-27d4-4d32-8021-80fc6b0d1472
# TODO
PlutoTest.@test ifft(arr_odd) ≈ ift(arr_odd)

# ╔═╡ a43b5bb5-d56b-43c7-aa96-ed2d9f217373
md"# 2 Convolution

From calculus class we know that a convolution can be expressed via Fourier transforms like

$U * h = \mathcal{F}^{-1}\left[\mathcal{F}[U] \cdot \mathcal{F}[h] \right] = \mathcal{F}^{-1}\left[\mathcal{F}[U] \cdot H \right]$

where $*$ is a convolution and $\cdot$ an elementwise multiplication. $H$ is the OTF of the PSF $h$.

Now implement it yourself!
"

# ╔═╡ a0362131-41ad-4151-b851-d75af22791e1
"""
	my_conv(U, h)

Calculates a FFT based convolution between U and h.
The output is a real array! 
So either use real valued transforms `rfft` or use `fft` with a final `real` conversion.

# Example
```julia-repl
julia> my_conv([1.0,0,0,0], [1.0, 0.5, 0.0, 0.5])
4-element Vector{ComplexF64}:
 1.0 + 0.0im
 0.5 + 0.0im
 0.0 + 0.0im
 0.5 + 0.0im
```
"""
function my_conv(U::AbstractArray{T, N}, h; center=ntuple(i -> 1, N)) where {T, N}
	# TODO
	return abs.(similar(U))
end

# ╔═╡ b568e681-7eaf-4c1a-9f36-f46163c35041
begin
	img = Float32.(testimage("mandril_gray"))
	Gray.(img)
end;

# ╔═╡ cdf36a7b-54ff-4ab7-9497-81d1ed684be3
# some simple kernel
kernel = IndexFunArrays.normal(img, sigma=7);

# ╔═╡ 0d9e2e69-e51f-483c-aadd-6bde1ab0c28a
md"You should see the blurry monkey here. If not, it might be wrong."

# ╔═╡ 5436a802-9d78-48ff-af79-2cf8065fd514
simshow(my_conv(img, kernel, center=(257, 257)))

# ╔═╡ 2eb0bea1-b6dc-417f-81b3-435fede2aa66
md"## 2 Test"

# ╔═╡ c7a4ce47-9bdf-46f3-88a4-888df19b9cb1
PlutoTest.@test my_conv([1.0,2,3,4], [1.0,0,0,0]) ≈ [1,2,3,4]

# ╔═╡ 3b479c44-4f6d-41ad-af70-04616a64154c
PlutoTest.@test my_conv([1.0,2,3,4], [0.0,1.0,0,0], center=2) ≈ [1,2,3,4]

# ╔═╡ 39807680-27cf-444d-9d93-d845baab61a4
PlutoTest.@test my_conv([1.0,2,3,4,5], [0.0,0.0,1.0,0,0], center=3) ≈ [1,2,3,4, 5]

# ╔═╡ 51bcf5d5-0405-467e-b528-275965fa1287
PlutoTest.@test my_conv(img, IndexFunArrays.delta(img, offset=(1,1))) ≈ img

# ╔═╡ da0d4ab1-b1d6-4817-ad38-6d3213c4e075
PlutoTest.@test my_conv([1, 2, 3, 4, 5, 6], [-1, 2, 2, 1, 3, 1]) ≈ [36.0, 32.0, 28.0, 30.0, 20.0, 22.0]

# ╔═╡ 09321d84-f9a7-412d-8fe4-3ece1bd90b21
md"# 3 Incoherent Image Formation

In this part we want to simulate an incoherent imaging situation (fluorescent microscope).
We simplify it by only considering an unitless `radius`.

As reminder, the qualitative procedure is as following:

We take a point source, Fourier transform it, take only the frequency inside the radius, go back to real space, take the absolute squared.

Pay attention where the frequencies should be located when you apply the `circ`.
"

# ╔═╡ ebfe7f3a-098a-41d8-a5bd-2c1e2b374fe0
begin
	"""
		circ(size, radius)
	
	`size` the size of the resulting array
	and `radius` the radius of the circle
	"""
	circ(size, radius) = rr2(size) .<= radius^2
end

# ╔═╡ 772fe3f9-6db8-4df4-8813-cbb334035351
"""
	calc_psf(arr_size, radius)

Calculate a simple PSF. `arr_size` is the output size and `radius` is the radius 
in Fourier space.
The output array is normalized by its sum to 1.
"""
function calc_psf(arr_size, radius)
	# TODO
	# return h
	return Array{Float32, 2}(undef, arr_size...) 
end

# ╔═╡ 33b3a409-2690-4a4a-a091-cdfe6a831c59
md"r = $(@bind r PlutoUI.Slider(0.01:0.1:10.0, show_value=true))"

# ╔═╡ f209c5de-c2bd-4a0b-bbfd-6831e0254023
simshow(calc_psf((64, 64), r))

# ╔═╡ 8b922c48-7d56-4e5b-b039-789e281c5fe1
md"r2 = $(@bind r2 PlutoUI.Slider(1:1:256, show_value=true))"

# ╔═╡ 119da3f0-0f1d-4a65-93cd-868f3c1d5f3e
h = calc_psf(size(img), r2);

# ╔═╡ 5830d782-67a3-4cae-847c-bdbdf0217aa7
# change this line such that the monkey is correctly centerd
# TODO
simshow(my_conv(img, h, center=size(h) .÷ 2 .+ 1))

# ╔═╡ fbdcb172-b079-46a5-b8f3-f5ece30fe25a
md"## 3 Test"

# ╔═╡ 66db7400-9791-4952-9812-39b22829b29a
# large radius is a perfect optical system -> delta peak
PlutoTest.@test calc_psf((2, 2), 1000) ≈  [0 0; 0 1]

# ╔═╡ baff7f36-d18a-4434-b979-662b6d44eb46
PlutoTest.@test sum(calc_psf((13, 12), 3)) ≈ 1

# ╔═╡ 7d92a5c5-d48d-429b-85fa-63904f21fc62
PlutoTest.@test minimum(calc_psf((13, 12), 3)) ≥ 0 

# ╔═╡ fd85fd4a-3e43-431c-aaa6-e92848c9e304
# compare to (approx) analytical solution
begin
	h2 = jinc.(7.25 / 9.219π * IndexFunArrays.rr((64, 64))).^2
	h2 ./= sum(h2)
	PlutoTest.@test ≈(1 .+ h2, 1 .+ calc_psf((64, 64), 7.25), rtol=0.001)
end

# ╔═╡ 08936039-7557-4ea3-8e26-5fbcdf12aec2
md"# 4 Generalized Wiener Filtering

A simpled deconvolution approach suited for Gaussian noise is the Wiener filter.
You can find the details in the slides.
Try to implement it here!
"

# ╔═╡ 7256f510-10ac-4fb2-99fc-0ffbcd45aae3
function wiener_filter(img, h, ϵ)
	# TODO
	return abs.(randn(size(img)))
end

# ╔═╡ 9f4d5c69-c5a8-4a34-9100-e8209ede71b4
begin
	# PSF
	h3 = ifftshift(calc_psf(size(img), 20))
	simshow(h3)
end

# ╔═╡ 2d5653cf-fa36-4a0e-99cd-4c9769b84705
begin
	img_b = my_conv(img, h3)
	simshow(img_b)
end

# ╔═╡ 0f8555c7-a1c2-4643-912b-6816019a848a
img_gauss = add_gauss(img_b, 0.1);

# ╔═╡ 34ced1ed-82f4-4c32-9cec-086ff2a62bce
img_poisson = poisson(img_b, 20);

# ╔═╡ cbd3d6c4-bf40-4431-a7b9-5b26b411a9b9
simshow([img_gauss img_poisson])

# ╔═╡ 29cc8e86-c109-411e-ba86-a8155f7c3a94
md"
pow1 = $(@bind pow1 Slider(-6:0.1:-0, show_value=true))

pow2 = $(@bind pow2 Slider(-6:0.1:-0, show_value=true))
"

# ╔═╡ 82542705-e24e-409e-a7bd-491ce750007e
img_gauss_wiener = wiener_filter(img_gauss, h3, 10^pow1);

# ╔═╡ e50febee-5c4b-44d7-9d78-29395a8c3ab6
img_poisson_wiener = wiener_filter(img_poisson, h3, 10^pow2);

# ╔═╡ 0b46da45-b743-40a1-a1f8-0581b6fd741a
simshow([img_gauss_wiener  img_poisson_wiener])

# ╔═╡ f85f964d-0443-4e19-b655-036f82a0ba69
md"## 4 Test"

# ╔═╡ 94f64888-b732-4367-a905-b57de684fcf7
PlutoTest.@test wiener_filter([1.0, 2.0], [1.0, 0.0], 0) ≈ [1.0, 2.0]

# ╔═╡ 3d89c0be-a888-4e2f-8820-85e07bd6be30
PlutoTest.@test  wiener_filter([1.0, 2.0], [1.0, 0.1], 0)  ≈ [0.808081, 1.91919] rtol=1e-5

# ╔═╡ c1e8419c-c748-490b-93d3-9bbcde4f4da9
md"# 5 Gradient Descent Optimization
In this part we want to implement an optimization routine with a _strange_ sensor.

The sensor has an additive Gaussian noise part but also an cubic gain behaviour.
See the function below
"

# ╔═╡ c337e8cf-bab9-46cd-aae8-22f6265f9cb7
begin
	function strange_sensor_f(value::T) where T
		value .^2
	end

	function strange_sensor(value::T) where T
		return abs.(randn(T) * 0.10f0 + strange_sensor_f(value))
	end
end

# ╔═╡ ff2288f0-2edf-4f22-90ef-af5777676ae7
# strange output
strange_sensor.([1.0 2; 3 4])

# ╔═╡ d23fe57a-7255-4d8f-a00a-eaec953213b2
md"
First we simulate the full `img` of the mandril with that sensor.
Clearly the appearance has changed due to the quadratic behaviour but also noise is visible.
"

# ╔═╡ a0f77be2-c40c-4372-a245-45b76ddc5861
begin
	img_strange = strange_sensor.(img) # todo
	simshow(img_strange)
end

# ╔═╡ cdc25e03-c6fd-4ee0-801e-24654a83d65d
md"## 5.1 Task
Since our sensor is quite noisy, we want to measure the image multiple times.
Try to complete `measure`.
"

# ╔═╡ 4af883c2-fdc9-458c-998b-ff0b8b5df146
"""
	measure(img, N)

Measure the `img` `N` times. Return a vector with the `N` measured images.
"""
function measure(img, N)
	# TODO
	return [img for i in 1:N]
end

# ╔═╡ c02969bc-3b5e-4422-a3c9-3a7d00cb3ddf
imgs = measure(img, 10); # TODO: simulate 10 times

# ╔═╡ e3a8ef81-0d19-4e75-b137-6454ce262991
simshow([reduce(hcat, imgs[1:end÷2]); reduce(hcat, imgs[end÷2+1:end])])

# ╔═╡ ebd74bcb-4187-469f-a206-ecf9041918f1
md"## 5.1 Test"

# ╔═╡ 7579c331-fa53-47f4-a479-312f2f7a3931
PlutoTest.@test typeof(measure([1.0], 12)) <: Vector

# ╔═╡ 3cca079c-1245-493b-ae6c-9db18e11835c
PlutoTest.@test length(measure([1.0], 12)) == 12

# ╔═╡ 0ae38ccc-9129-4c61-9810-dbc448e8a07f
begin
	Random.seed!(42)
	a = measure([1.0 2.0; 3 4], 3)
end;

# ╔═╡ 4f5a9a27-1893-49e4-a2e9-59ba3a2340a2
begin
		Random.seed!(42)
		b = [abs.([1.0 2.0; 3 4].^2 .+ randn((2,2)) * 0.1f0) for i = 1:3]
end;

# ╔═╡ 37ce2f96-377d-48d7-bc56-f83b0cce349c
PlutoTest.@test a ≈ b

# ╔═╡ cd0829a1-cdf1-45c9-a4dd-91bb3ec5bb03
md"## 5.2 Loss Function
Having our 10 images we would like to retrieve a best guess for the underlying image.
Taking the mean does not work in this case since the input image is modified by `strange_sensor_f`.

Therefore, we interprete the reconstruction as an optimization problem.
Try to find the corresponding pages in the lecture.

Generally, we want to minimize a loss function which looks like

$\underset{\mu}{\mathrm{argmin}}\, \mathcal L = \sum_{\text{img} \, \in \, \text{imgs}} \sum_{p_i \, \in \, \text{img}} (\text{img}[i] - f(\mu[i]))^2$


So we sum over all measured images. For each image, we additionally sum each pixel and  compare it with $f(\mu[i])$.
$f$ is the same as `strange_sensor_f` and by calling $f(\mu[i])$ we apply $f$ to the reconstruction $\mu$. Via that function call, we hope that we find a $\mu$ which fits to the measurments. So $f$ is the forward model of the sensor (without the noise part).
"

# ╔═╡ 328b5994-b1c4-4850-b8d3-2f3781bed99f
"""
	loss(imgs, μ)

Calculate the loss between the `ìmgs` and `μ`.
Basically implement the sum of the square difference value.
Don't forget to apply `strange_sensor_f` to `μ` before!

Using two for loops is perfectly fine!
"""
function loss(imgs::Vector{<:AbstractArray{T, N}}, μ::AbstractArray{T, N}) where {T, N}
	# TODO
	return rand()
end

# ╔═╡ 78da1c47-1de7-4d95-a1e1-cd60e912a0fc
# comparison with ground truth image
loss(imgs, img) 

# ╔═╡ 3d182fa9-84bb-4504-8945-bf653ce1f99d
# the ground truth is not the minimum...
# why not?
loss(imgs, img .+0.0001f0) 

# ╔═╡ b56a311e-912a-4e7c-b821-4124b847194a
md"## 5.2 Test"

# ╔═╡ 840ddb7e-56ba-4d7a-8f08-cee67128315f
PlutoTest.@test loss([[1]], [2]) isa Number

# ╔═╡ 528194f2-9f74-4dda-94f8-3e4f0e515bc9
PlutoTest.@test loss([[1]], [2])  ≈ 9

# ╔═╡ a203fb4f-c3d6-4d46-9c58-2454b71b0e1b
PlutoTest.@test loss([[1], [2], [4.0]], sqrt.([2]))  ≈ 5

# ╔═╡ 063ef53f-d152-4282-af5f-7f2addce3ab0
md"## 5.3 Gradient of `strange_sensor_f`
For the optimization we later want to apply a gradient descent optimization scheme.
Hence, we need a few gradients!
"

# ╔═╡ 95392d70-07cf-47a7-8707-f1dab177e7a5
"""
	 grad_strange_sensor_f(value::T)

Calculate the gradient of `strange_sensor_f`.
Use the rules you know already from school!
"""
function grad_strange_sensor_f(value::T) where T
	# TODO
	return rand()
end

# ╔═╡ 01de98fe-5537-4745-acf8-99daa0b8aa6f
grad_strange_sensor_f(1)

# ╔═╡ 75549871-14a3-4b9f-9b98-cccc34ed315d
md"## 5.3 Test"

# ╔═╡ 867b3bf0-1208-48a1-b1ee-904d11b28e1f
PlutoTest.@test grad_strange_sensor_f(0) ≈ 0

# ╔═╡ 7c2d2d6d-f933-452e-9870-7dce7ad4bf1d
# comparison with automatic differentation package!
PlutoTest.@test ForwardDiff.derivative(strange_sensor_f, 42) ≈ grad_strange_sensor_f(42)

# ╔═╡ cd363a92-a473-47f8-b41b-c5d5249fec90
md"## 5.4 Gradient of Loss
Now we come to the last part of the gradient.
We need the full gradient of the loss function.

The gradient of the loss function with respect to `μ` will be an array again.
That is plausible since the loss function accounts for all pixels. By changing a single pixel we also change the value of the loss function. Therefore, for each pixel a gradient exists describing the influence to the loss value.

You need to derive the loss with respect to $\mu$. But don't forget about the chain rule for `strange_sensor_f`!

Again, two for loops are perfectly fine!
"

# ╔═╡ 1ad4cf95-a3fb-4e45-b7bd-80ab724c158e
function gradient(imgs::Vector{<:AbstractArray{T, N}}, μ::AbstractArray{T, N}) where {T, N}
	# TODO
	return similar(imgs[1])
end

# ╔═╡ cb41a2f9-a0c5-4eeb-af51-9bf40a5d7fd6
gradient(imgs, img)

# ╔═╡ f72450cb-a54e-4827-9c33-7a44f561bd43
md"## 5.4 Test"

# ╔═╡ 7e6eb33d-0721-4bb6-a9ee-1cfb07a17e46
PlutoTest.@test gradient([[1.0], [2.0]], [2.5]) ≈ [95]

# ╔═╡ cf550002-13ce-4788-8206-790fb355a91b
PlutoTest.@test gradient([[3.5, 2.1], [2.0, 0.0]], [3.23, 23.2]) isa Vector

# ╔═╡ bfd35c9f-1449-4d86-bb50-9bbbd6ea2c30
PlutoTest.@test size(gradient([[3.5, 2.1], [2.0, 0.0]], [3.23, 23.2])) == (2,)

# ╔═╡ c1795dd9-f8a5-472d-bd82-7871e8603534
PlutoTest.@test  gradient([[3.5, 2.1], [2.0, 0.0]], [3.23, 23.2]) ≈ [198.526136, 99702.46399999999]

# ╔═╡ c5436e99-513b-4c51-b3ad-f82318304a3e
md"## 5.5 Gradient Descent"

# ╔═╡ b393a0f5-af03-4ded-8e22-886022f5da30
"""
	gradient_descent(imgs, μ, N_iter, step_size)

Runs a gradient descent using `gradient` and the experimental images `imgs`.
`N_iter` is the number of iterations, `step_size` is the step size of the gradient step.

The optimized output is `μ_optimized`.
"""
function gradient_descent(imgs, μ, N_iter, step_size)
	# TODO
	μ_optimized = copy(μ)
	for i = 1:N_iter
		μ_optimized .= μ_optimized .- step_size .* gradient(imgs, μ_optimized)
	end
	return μ_optimized
end

# ╔═╡ c81c5eaa-0e0f-4730-9148-ec0e7b63cdd4
μ_init = similar(imgs[1]) # TODO, what could be a good initilization?

# ╔═╡ 32f6d83a-d5db-4a21-b2c4-d8876de38c46
md"
Try to change the gradient step and the number of iterations.

Number of iterations $(@bind N_iter Slider(1:100, show_value=true))

step size $(@bind step_size Slider(1f-5:1f-5:1f-1, show_value=true))
"

# ╔═╡ f7d7cf86-ece9-406b-af78-f86177a767c4
μ_optimized = gradient_descent(imgs, μ_init, N_iter, step_size);

# ╔═╡ 2110469f-16f9-4454-996a-a949a94fffa3
# that value should get smaller with more iterations
# smaller is better
loss(imgs, μ_optimized)

# ╔═╡ 3643e55e-3205-4561-b271-058d59f46685
# this value should be larger than ↑
loss(imgs, μ_init)

# ╔═╡ 85c81f80-c959-4682-932b-26fe7f415f4d
Gray.(μ_optimized)

# ╔═╡ 09313c1d-102f-4e0b-b9c8-7d2e904a1dd6
md"## 5.5 Test"

# ╔═╡ b2e7aae8-2777-4657-ba3d-f0dfcb581ede
PlutoTest.@test gradient_descent([[1.0], [2.0]], [3.0], 3, 0.01) ≈ [1.21021] rtol=1f-4

# ╔═╡ ae027113-998b-467e-8967-cc416d480377
PlutoTest.@test loss(imgs, μ_init) > loss(imgs, μ_optimized)

# ╔═╡ 4ff0afae-af4c-45d0-83b2-4fb5eead2fa1
PlutoTest.@test loss(imgs, gradient_descent(imgs, (imgs[1] .+ imgs[2]) / 2, 15, 0.000001);) < loss(imgs, imgs[1])

# ╔═╡ 2374bfc1-cd14-4bee-b5e5-bb15e0e6e253
begin
	μ_mean = sqrt.(reduce((a,b) -> a .+ b, imgs) ./ 10)
	Gray.(μ_mean)
end;

# ╔═╡ 9360c9a0-6a0f-4ff6-a215-23ce5c3c8099
md"#### Can you beat the _mean_?"

# ╔═╡ a837b22c-61f1-4a98-a362-c08b405c4dca
PlutoTest.@test loss(imgs, μ_optimized) < loss(imgs, μ_mean)

# ╔═╡ e2aec734-ce2b-4761-b3d0-b559df8b17da
md"Note that in this particular example the mean can be proven to be the best estimator, so it is not surprising to not be able to beat the mean. Yet, if the model cannot be inverted, we still may want to use a gradient-based optimization to solve the problem."

# ╔═╡ 959b4745-7742-459a-bef0-e374ec4aec17
[simshow(img) simshow(μ_mean) simshow(μ_optimized)]

# ╔═╡ b3d54170-466c-4b0d-a337-280ef4ea87f3
md"## 5.6 Test
To check wether your gradient and your loss is correct, see also the following output of Optim (a sophisticated package for optimization).
If the output is not a nice Mandrill, there is most likely something wrong with your loss/gradients.
"

# ╔═╡ eafe5f6e-18ab-4341-8435-a3c0b9423b35
begin
	g_optim!(G, x) = isnothing(G) ? nothing : G .= gradient(imgs, x)
	f_optim(x) = loss(imgs, x)
end

# ╔═╡ 13304ea6-0639-4863-827e-d66a8630bc63
begin
	μ_optim_init = copy(μ_init)
	res = optimize(f_optim, g_optim!, μ_optim_init, ConjugateGradient(), Optim.Options(iterations=50))
	μ_optim = Optim.minimizer(res)
	res
end

# ╔═╡ c685ac28-79e8-4eda-8c0f-0944a1276691
simshow(Optim.minimizer(res))

# ╔═╡ 7a691c80-e08d-423a-a69a-572f02fdddda
loss(imgs, μ_optim)

# ╔═╡ 19a0e512-46d6-429e-93f0-05e2faf95218
md"## Can you beat Optim?"

# ╔═╡ abe5aca2-61ce-4f99-ad75-0c84293cd7b3
PlutoTest.@test loss(imgs, μ_optimized) < loss(imgs, μ_optim)

# ╔═╡ 1171004b-96a1-4ee4-bbd5-2b7f2b1d582a
# yes, but hard with an hand optimized iterations/value pair

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
FourierTools = "b18b359b-aebc-45ac-a139-9c0ccbb2871e"
ImageIO = "82e4d734-157c-48bb-816b-45c225c6df19"
ImageShow = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
IndexFunArrays = "613c443e-d742-454e-bfc6-1d7f8dd76566"
Noise = "81d43f40-5267-43b7-ae1c-8b967f377efa"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
TestImages = "5e47fb64-e119-507b-a336-dd2b206d9990"

[compat]
FFTW = "~1.5.0"
ForwardDiff = "~0.10.33"
FourierTools = "~0.3.7"
ImageIO = "~0.6.6"
ImageShow = "~0.3.6"
IndexFunArrays = "~0.2.5"
Noise = "~0.3.2"
Optim = "~1.7.4"
Plots = "~1.36.6"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.49"
SpecialFunctions = "~2.1.7"
TestImages = "~1.7.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.3"
manifest_format = "2.0"
project_hash = "eed75f377c1c0cccf6b70ac1c218877c54cf8509"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[deps.AbstractNFFTs]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "292e21e99dedb8621c15f185b8fdb4260bb3c429"
uuid = "7f219486-4aa7-41d6-80a7-e08ef20ceed7"
version = "0.8.2"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "c46fb7dd1d8ca1d213ba25848a5ec4e47a1a1b08"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.26"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "1dd4d9f5beebac0c03446918741b1a03dc5e5788"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.6"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "7fe6d92c4f281cf4ca6f2fba0ce7b299742da7ca"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.37"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.BasicInterpolators]]
deps = ["LinearAlgebra", "Memoize", "Random"]
git-tree-sha1 = "3f7be532673fc4a22825e7884e9e0e876236b12a"
uuid = "26cce99e-4866-4b6d-ab74-862489e035e0"
version = "0.7.1"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "00a2cccc7f098ff3b66806862d275ca3db9e6e5a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.ContextVariablesX]]
deps = ["Compat", "Logging", "UUIDs"]
git-tree-sha1 = "25cc3803f1030ab855e383129dcd3dc294e322cc"
uuid = "6add18c4-b38d-439d-96f6-d6bc489c04c5"
version = "0.1.3"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "e08915633fcb3ea83bf9d6126292e5bc5c739922"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.13.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "c5b6685d53f933c11404a3ae9822afe30d522494"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.12.2"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "c36550cb29cbe373e95b3f40486b9a4148f89ffd"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.2"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FLoops]]
deps = ["BangBang", "Compat", "FLoopsBase", "InitialValues", "JuliaVariables", "MLStyle", "Serialization", "Setfield", "Transducers"]
git-tree-sha1 = "ffb97765602e3cbe59a0589d237bf07f245a8576"
uuid = "cc61a311-1640-44b5-9fba-1b764f453329"
version = "0.2.1"

[[deps.FLoopsBase]]
deps = ["ContextVariablesX"]
git-tree-sha1 = "656f7a6859be8673bf1f35da5670246b923964f7"
uuid = "b9860ae5-e623-471e-878b-f6a53c775ea6"
version = "0.1.1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "802bfc139833d2ba893dd9e62ba1767c88d708ae"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.5"

[[deps.FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "04ed1f0029b6b3af88343e439b995141cb0d0b8d"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.17.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "10fa12fe96e4d76acfa738f4df2126589a67374f"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.33"

[[deps.FourierTools]]
deps = ["ChainRulesCore", "FFTW", "IndexFunArrays", "LinearAlgebra", "NDTools", "NFFT", "PaddedViews", "Reexport", "ShiftedArrays"]
git-tree-sha1 = "3cbae3e75b991e2b58034af6de389a6d26e33cfe"
uuid = "b18b359b-aebc-45ac-a139-9c0ccbb2871e"
version = "0.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "051072ff2accc6e0e87b708ddee39b18aa04a0bc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.71.1"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "501a4bf76fd679e7fcd678725d5072177392e756"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.71.1+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78e2c69783c9753a91cdae88a8d432be85a2ab5e"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fb83fbe02fe57f2c068013aa94bcdf6760d3a7a7"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+1"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "e1acc37ed078d99a714ed8376446f92a5535ca65"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.5.5"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "c54b581a83008dc7f292e205f4c409ab5caa0f04"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.10"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "342f789fd041a55166764c351da1710db97ce0e0"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.6"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils", "Libdl", "Pkg", "Random"]
git-tree-sha1 = "5bc1cb62e0c5f1005868358db0692c994c3a13c6"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.2.1"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "124626988534986113cfd876e3093e4a03890f58"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.12+3"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "36cbaebed194b292590cba2593da27b34763804a"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.8"

[[deps.ImageShow]]
deps = ["Base64", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "b563cf9ae75a635592fc73d3eb78b86220e55bd8"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.6"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.IndexFunArrays]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "c7e4b47fa1cd2761794b96b3e6ac1d7a0c2133aa"
uuid = "613c443e-d742-454e-bfc6-1d7f8dd76566"
version = "0.2.5"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "16c0cc91853084cb5f58a78bd209513900206ce6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.JuliaVariables]]
deps = ["MLStyle", "NameResolution"]
git-tree-sha1 = "49fb3cb53362ddadb4415e9b73926d6b40709e70"
uuid = "b14d175d-62b4-44ba-8fb7-3064adc8c3ec"
version = "0.2.4"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "946607f84feb96220f480e0422d3484c49c00239"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.19"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MLStyle]]
git-tree-sha1 = "060ef7956fef2dc06b0e63b294f7dbfbcbdc7ea2"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.16"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Memoize]]
deps = ["MacroTools"]
git-tree-sha1 = "2b1dfcba103de714d31c033b5dacc2e4a12c7caa"
uuid = "c03570c3-d221-55d1-a50c-7939bbd78826"
version = "0.4.4"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "4d5917a26ca33c66c8e5ca3247bd163624d35493"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NDTools]]
deps = ["LinearAlgebra", "OffsetArrays", "PaddedViews", "Random", "Statistics"]
git-tree-sha1 = "f9de89bade7fa11d65360133b6ea0cbe917935e7"
uuid = "98581153-e998-4eef-8d0d-5ec2c052313d"
version = "0.5.1"

[[deps.NFFT]]
deps = ["AbstractNFFTs", "BasicInterpolators", "Distributed", "FFTW", "FLoops", "LinearAlgebra", "Printf", "Random", "Reexport", "SnoopPrecompile", "SparseArrays", "SpecialFunctions"]
git-tree-sha1 = "93a5f32dd6cf09456b0b81afcb8fc29f06535ffd"
uuid = "efe261a4-0d2b-5849-be55-fc731d526b0d"
version = "0.13.3"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NameResolution]]
deps = ["PrettyPrint"]
git-tree-sha1 = "1a0fa0e9613f46c9b8c11eee38ebb4f590013c5e"
uuid = "71a1bf82-56d0-4bbc-8a3c-48b961074391"
version = "0.1.5"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "5ae7ca23e13855b3aba94550f26146c01d259267"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Noise]]
deps = ["ImageCore", "PoissonRandom", "Random"]
git-tree-sha1 = "1427315f223bc7c754c1d97a12c2b5fc059dafbc"
uuid = "81d43f40-5267-43b7-ae1c-8b967f377efa"
version = "0.3.2"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "f71d8950b724e9ff6110fc948dff5a329f901d64"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "df6830e37943c7aaa10023471ca47fb3065cc3c4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "1903afc76b7d01719d9c30d3c7d501b61db96721"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.4"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "f809158b27eba0c18c269cf2a2be6ed751d3e81d"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.17"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "b64719e8b4504983c7fca6cc9db3ebc8acc2a4d6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "21303256d239f6b484977314674aef4bb1fe4420"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "6a9521b955b816aa500462951aa67f3e4467248a"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.36.6"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "45f9da1ceee5078267eb273d065e8aa2f2515790"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.3"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyPrint]]
git-tree-sha1 = "632eb4abab3449ab30c5e1afaa874f0b98b586e4"
uuid = "8162dcfd-2161-5ef2-ae6c-7681170c5f98"
version = "0.2.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "18c35ed630d7229c5584b945641a73ca83fb5213"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.2"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "e974477be88cb5e3040009f3767611bc6357846f"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.11"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.ShiftedArrays]]
git-tree-sha1 = "503688b59397b3307443af35cd953a13e8005c16"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "2.0.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "ffc098086f35909741f71ce21d03dadf0d2bfa76"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.11"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StringDistances]]
deps = ["Distances", "StatsAPI"]
git-tree-sha1 = "ceeef74797d961aee825aabf71446d6aba898acb"
uuid = "88034a9c-02f8-509d-84a9-84ec65e18404"
version = "0.11.2"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TestImages]]
deps = ["AxisArrays", "ColorTypes", "FileIO", "ImageIO", "ImageMagick", "OffsetArrays", "Pkg", "StringDistances"]
git-tree-sha1 = "03492434a1bdde3026288939fc31b5660407b624"
uuid = "5e47fb64-e119-507b-a336-dd2b206d9990"
version = "1.7.1"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "f8cd5b95aae14d3d88da725414bdde342457366f"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.2"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "8a75929dcd3c38611db2f8d08546decb514fcadf"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.9"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "c42fa452a60f022e9e087823b47e5a5f8adc53d5"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.75"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╠═db61ecb1-6ed9-44d9-b8e4-b50622efbac1
# ╠═fbc3ec99-0c73-40e7-bd21-933bb5a52755
# ╟─82726a23-060c-4c7d-91c1-e78a22fb98d6
# ╠═ebe90686-c72f-4a42-a644-f06babe00f41
# ╠═5053c100-5f27-4f2f-acbe-e12c9070feac
# ╠═361a7a34-d70b-4597-948c-80ed913a562b
# ╟─653a94c9-cb89-47df-8ba5-0ce683440b78
# ╠═b65e0745-82c6-48e1-a956-fa71583dad15
# ╟─c53325bd-eeb6-4663-867d-44fd3ace273b
# ╠═b57917e8-3d2e-433b-84df-6129f760f954
# ╟─6e0e223d-9107-4121-880a-8ed5d7e5e9c2
# ╠═5213f82b-5e06-4cee-a0d1-21f1c5cc2998
# ╠═64103e19-ce79-48b7-bac9-b429f6d423c2
# ╟─f380db27-8773-413b-978f-4496fc585ae3
# ╠═0464d972-fcaf-4f83-be7d-073d217e8f4c
# ╠═87fabbfe-e25b-467a-9fbd-ceec98ba4ed6
# ╠═c537e21b-5ed2-4cac-947a-e052474a2442
# ╠═3093ae8d-27d4-4d32-8021-80fc6b0d1472
# ╟─a43b5bb5-d56b-43c7-aa96-ed2d9f217373
# ╠═a0362131-41ad-4151-b851-d75af22791e1
# ╠═b568e681-7eaf-4c1a-9f36-f46163c35041
# ╠═cdf36a7b-54ff-4ab7-9497-81d1ed684be3
# ╟─0d9e2e69-e51f-483c-aadd-6bde1ab0c28a
# ╠═5436a802-9d78-48ff-af79-2cf8065fd514
# ╟─2eb0bea1-b6dc-417f-81b3-435fede2aa66
# ╠═c7a4ce47-9bdf-46f3-88a4-888df19b9cb1
# ╠═3b479c44-4f6d-41ad-af70-04616a64154c
# ╠═39807680-27cf-444d-9d93-d845baab61a4
# ╠═51bcf5d5-0405-467e-b528-275965fa1287
# ╠═da0d4ab1-b1d6-4817-ad38-6d3213c4e075
# ╟─09321d84-f9a7-412d-8fe4-3ece1bd90b21
# ╠═ebfe7f3a-098a-41d8-a5bd-2c1e2b374fe0
# ╠═772fe3f9-6db8-4df4-8813-cbb334035351
# ╟─33b3a409-2690-4a4a-a091-cdfe6a831c59
# ╠═f209c5de-c2bd-4a0b-bbfd-6831e0254023
# ╟─8b922c48-7d56-4e5b-b039-789e281c5fe1
# ╠═119da3f0-0f1d-4a65-93cd-868f3c1d5f3e
# ╠═5830d782-67a3-4cae-847c-bdbdf0217aa7
# ╟─fbdcb172-b079-46a5-b8f3-f5ece30fe25a
# ╠═66db7400-9791-4952-9812-39b22829b29a
# ╠═baff7f36-d18a-4434-b979-662b6d44eb46
# ╠═7d92a5c5-d48d-429b-85fa-63904f21fc62
# ╠═fd85fd4a-3e43-431c-aaa6-e92848c9e304
# ╟─08936039-7557-4ea3-8e26-5fbcdf12aec2
# ╠═7256f510-10ac-4fb2-99fc-0ffbcd45aae3
# ╠═9f4d5c69-c5a8-4a34-9100-e8209ede71b4
# ╠═2d5653cf-fa36-4a0e-99cd-4c9769b84705
# ╠═0f8555c7-a1c2-4643-912b-6816019a848a
# ╠═34ced1ed-82f4-4c32-9cec-086ff2a62bce
# ╠═cbd3d6c4-bf40-4431-a7b9-5b26b411a9b9
# ╠═29cc8e86-c109-411e-ba86-a8155f7c3a94
# ╠═82542705-e24e-409e-a7bd-491ce750007e
# ╠═e50febee-5c4b-44d7-9d78-29395a8c3ab6
# ╠═0b46da45-b743-40a1-a1f8-0581b6fd741a
# ╟─f85f964d-0443-4e19-b655-036f82a0ba69
# ╠═94f64888-b732-4367-a905-b57de684fcf7
# ╠═3d89c0be-a888-4e2f-8820-85e07bd6be30
# ╟─c1e8419c-c748-490b-93d3-9bbcde4f4da9
# ╠═c337e8cf-bab9-46cd-aae8-22f6265f9cb7
# ╠═ff2288f0-2edf-4f22-90ef-af5777676ae7
# ╟─d23fe57a-7255-4d8f-a00a-eaec953213b2
# ╠═a0f77be2-c40c-4372-a245-45b76ddc5861
# ╟─cdc25e03-c6fd-4ee0-801e-24654a83d65d
# ╠═4af883c2-fdc9-458c-998b-ff0b8b5df146
# ╠═c02969bc-3b5e-4422-a3c9-3a7d00cb3ddf
# ╠═e3a8ef81-0d19-4e75-b137-6454ce262991
# ╟─ebd74bcb-4187-469f-a206-ecf9041918f1
# ╠═7579c331-fa53-47f4-a479-312f2f7a3931
# ╠═3cca079c-1245-493b-ae6c-9db18e11835c
# ╟─0ae38ccc-9129-4c61-9810-dbc448e8a07f
# ╟─4f5a9a27-1893-49e4-a2e9-59ba3a2340a2
# ╠═37ce2f96-377d-48d7-bc56-f83b0cce349c
# ╟─cd0829a1-cdf1-45c9-a4dd-91bb3ec5bb03
# ╠═328b5994-b1c4-4850-b8d3-2f3781bed99f
# ╠═78da1c47-1de7-4d95-a1e1-cd60e912a0fc
# ╠═3d182fa9-84bb-4504-8945-bf653ce1f99d
# ╟─b56a311e-912a-4e7c-b821-4124b847194a
# ╠═840ddb7e-56ba-4d7a-8f08-cee67128315f
# ╠═528194f2-9f74-4dda-94f8-3e4f0e515bc9
# ╠═a203fb4f-c3d6-4d46-9c58-2454b71b0e1b
# ╟─063ef53f-d152-4282-af5f-7f2addce3ab0
# ╠═95392d70-07cf-47a7-8707-f1dab177e7a5
# ╠═01de98fe-5537-4745-acf8-99daa0b8aa6f
# ╟─75549871-14a3-4b9f-9b98-cccc34ed315d
# ╠═867b3bf0-1208-48a1-b1ee-904d11b28e1f
# ╠═7c2d2d6d-f933-452e-9870-7dce7ad4bf1d
# ╟─cd363a92-a473-47f8-b41b-c5d5249fec90
# ╠═1ad4cf95-a3fb-4e45-b7bd-80ab724c158e
# ╠═cb41a2f9-a0c5-4eeb-af51-9bf40a5d7fd6
# ╟─f72450cb-a54e-4827-9c33-7a44f561bd43
# ╠═7e6eb33d-0721-4bb6-a9ee-1cfb07a17e46
# ╠═cf550002-13ce-4788-8206-790fb355a91b
# ╠═bfd35c9f-1449-4d86-bb50-9bbbd6ea2c30
# ╠═c1795dd9-f8a5-472d-bd82-7871e8603534
# ╟─c5436e99-513b-4c51-b3ad-f82318304a3e
# ╠═b393a0f5-af03-4ded-8e22-886022f5da30
# ╠═c81c5eaa-0e0f-4730-9148-ec0e7b63cdd4
# ╟─32f6d83a-d5db-4a21-b2c4-d8876de38c46
# ╠═f7d7cf86-ece9-406b-af78-f86177a767c4
# ╠═2110469f-16f9-4454-996a-a949a94fffa3
# ╠═3643e55e-3205-4561-b271-058d59f46685
# ╠═85c81f80-c959-4682-932b-26fe7f415f4d
# ╟─09313c1d-102f-4e0b-b9c8-7d2e904a1dd6
# ╠═b2e7aae8-2777-4657-ba3d-f0dfcb581ede
# ╠═ae027113-998b-467e-8967-cc416d480377
# ╠═4ff0afae-af4c-45d0-83b2-4fb5eead2fa1
# ╠═2374bfc1-cd14-4bee-b5e5-bb15e0e6e253
# ╟─9360c9a0-6a0f-4ff6-a215-23ce5c3c8099
# ╠═a837b22c-61f1-4a98-a362-c08b405c4dca
# ╟─e2aec734-ce2b-4761-b3d0-b559df8b17da
# ╠═959b4745-7742-459a-bef0-e374ec4aec17
# ╟─b3d54170-466c-4b0d-a337-280ef4ea87f3
# ╠═eafe5f6e-18ab-4341-8435-a3c0b9423b35
# ╠═13304ea6-0639-4863-827e-d66a8630bc63
# ╠═c685ac28-79e8-4eda-8c0f-0944a1276691
# ╠═7a691c80-e08d-423a-a69a-572f02fdddda
# ╟─19a0e512-46d6-429e-93f0-05e2faf95218
# ╠═abe5aca2-61ce-4f99-ad75-0c84293cd7b3
# ╟─1171004b-96a1-4ee4-bbd5-2b7f2b1d582a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
