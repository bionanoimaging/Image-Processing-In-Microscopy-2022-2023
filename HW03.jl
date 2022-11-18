### A Pluto.jl notebook ###
# v0.19.14

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

# ╔═╡ 77e6389e-3491-45b6-9f76-459b15a8922d
using TestImages, ImageShow, Random, PlutoTest, Plots, FourierTools, Statistics, Colors

# ╔═╡ ee4b6b2d-119e-4ac3-a1f5-500ff045bfc5
using FFTW, PlutoUI

# ╔═╡ c7150c2c-afeb-4158-b381-3341bb678fd5
using Noise

# ╔═╡ 31c31237-c7bd-4fdb-8f7f-219d323f80e7
"""
    my_show(arr::AbstractArray{<:Real}; set_one=false, set_zero=false)
Displays a real valued array . Brightness encodes magnitude.
Works within Jupyter and Pluto.
## Keyword args
* `set_one=false` divides by the maximum to set maximum to 1
* `set_zero=false` subtracts the minimum to set minimum to 1
"""
function gray_show(arr::AbstractArray{<:Real}; set_one=true, set_zero=false)
    arr = set_zero ? arr .- minimum(arr) : arr
    arr = set_one ? arr ./ maximum(arr) : arr
    Gray.(arr)
end

# ╔═╡ 7148e705-8b90-46ea-8763-778897daf9a7
md"# Homework 03

In this homework we cover several topics like image correction, Fourier transforms, Fourier Shifting and Sampling.
"

# ╔═╡ 83d9c252-2a8d-47a3-b247-a62b6a175201
md"## 1. Flat Fielding
As discussed in the lecture, raw images might be affected by dirty optics or sensor errors.
In this exercise we want to apply such a correction to simulated images.
"

# ╔═╡ 39c623c7-6063-4360-a2b7-95ea42c5262d
begin
	img = Float32.(testimage("mandril_gray.tif"))
	gray_show(img)
end

# ╔═╡ eacd8fdf-971c-46ca-914f-f73f65e28ca8
md"### Simulate Deteriorated Image
In the first step, we degrade our image with some artifacts.
We add an offset to our image which is caused by sensor artifacts.
Second, we add some dirty spots to the image which are for example caused by dust on the sensor.
"

# ╔═╡ 25dde5cb-9ed6-498e-b511-50ed0640d0c0
"""
	simulate_bad_microscope(img)

The ideal `img` is degraded by artifacts.
"""
function simulate_bad_microscope(img::AbstractArray{T}) where T
	# generate same random numbers!
	rng = MersenneTwister(2021);
	# some weird dark field slope offset
	dark = zero.(img)
	dark .+= range(0,0.3, length=size(img, 1))

	xs = 1:size(img, 1)

	ff = one.(img)
	for i = 1:20 # emulates some dust spots as modulations to the flatfield ff
		x, y = rand(rng, xs), rand(rng, xs)
		sca = abs(randn(rng) .* T(0.2))
		r = sqrt.((xs .- x).^2 .+ (xs' .- y).^2)
		ff .-= T(0.2) .* sinc.(r .* sca)
	end

	return abs.(dark .+ ff .* img)
end

# ╔═╡ a6c024a4-7cb2-4f8d-be24-8727d75679da
img_dirty = simulate_bad_microscope(img);

# ╔═╡ d5c127ee-1ecf-4ed0-9cd9-6fcdefa95847
gray_show(img_dirty)

# ╔═╡ 6e95cbcb-0d89-4eca-8599-b266076708a7
md"## 1.1 Task 
Extract both a flat field and a dark field image with `simulate_bad_microscope`.
Use the function `simulate_bad_microscope` to obtain a measured flat field and measured dark field.
"

# ╔═╡ 910cb7fa-b1fc-4106-acc4-3bf56edfb83d
flat_field = similar(img); # TODO

# ╔═╡ cbebb311-c521-4fa6-a909-1a4df8da4d2d
gray_show(flat_field)

# ╔═╡ 8bb737cb-2d59-49ba-9b89-afef6ca7bbf8
dark_field = similar(img); # TODO

# ╔═╡ 52bc8b3d-542e-4a35-b38c-35840f4075d4
gray_show(dark_field)

# ╔═╡ e3448ffe-5fbf-4bed-bf03-6eb9031e6c12
"""
	correct_image(img_dirty, dark_field, flat_field)

Correct a dirty image `img_dirty` by providing the `dark_field` and `flat_field`.

Conceptually we divide by the `flat_field` since the measured `flat_field`
describes the deviations from an ideal white image.
The `dark_field` has to be subtracted before since that only contains sensor artifacts.
"""
function correct_image(img_dirty, dark_field, flat_field)
	# TODO
	return img_dirty
end

# ╔═╡ d9308666-efaa-4fcc-a083-f3d12dbab1e2
img_fixed = similar(img_dirty); # todo

# ╔═╡ 493b6bb4-1dfc-4a82-9360-40944b29826c
gray_show(img_fixed)

# ╔═╡ c80a2ce8-6620-4de1-bb78-5af97caf466d
md"### 1.1 Test"

# ╔═╡ e2e470e4-ff73-4dd5-973c-14e1535b01ea
PlutoTest.@test img_fixed ≈ img

# ╔═╡ 0a957f4a-a975-4087-b458-3236940eda43
PlutoTest.@test img ≈ correct_image(img, zero.(img), one.(img))

# ╔═╡ af9866f4-2b03-4dc4-9f8a-e02ae46e2f9c
md"# 2 Sampling
If a signal is not properly sampled, we can see a lot of artifacts which are called aliasing.
Try to determine experimentally (via the slider) which is the correct value for Nyquist sampling.

Afterwards, calculate the correct sampling in the variable $N$!

Don't state the pure number but rather a short equation which Julia calculates (something like `1 * 2 * 3`).
Do you see a difference between your expected and your calculated value?
"

# ╔═╡ 7a481af3-4ec8-4a58-bb84-e9d3d5381570
md"N=$(@bind N Slider(2:300, show_value=true))"

# ╔═╡ 4006ccf7-df8e-440c-8b3c-af525bb50ed4
begin
	xs = range(0, 8, length=N+1)[begin:end-1]
	f(x) =  sin(x * 2π * 5)
end

# ╔═╡ c871cdbd-b560-46f2-918c-45d24d2f9b75
plot(xs, f.(xs))

# ╔═╡ f6b8cd92-c1ab-445b-8b4f-04f521b72cad
N_correct = 2 # TODO

# ╔═╡ 51b2d7be-fe94-44e9-8e67-73aafb636ef1
begin
	xs_correct = range(0, 8, length=round(Int, N_correct + 1))[begin:end-1]
	plot(xs_correct, f.(xs_correct))
end

# ╔═╡ a19c8850-9a7d-4e26-908b-a067a9006e84
md"No tests, since they would reveal the answer." 

# ╔═╡ b0c67bd6-b5ef-4359-8727-1a3bfefcf759
md"# 3 FFTs

Another very important topics are FFTs!
The standard library in many languages is based on FFTW.
In Julia FFTW.jl is the recommended library. 

The discrete Fourier transform is given by

$$X_k = \sum_{n=0}^{N-1} x_n \exp\left(-i2 \pi \frac{kn}{N} \right)$$
and the inverse by

$$x_n = \frac1{N}\sum_{k=0}^{N-1} X_k \exp\left(i2 \pi \frac{kn}{N} \right)$$.

The FFT algorithm evaluates those equations very efficiently ($\mathcal O(N \cdot \log N)$) by avoiding the explicit sum ($\mathcal O(N^2)$).

Further, the FFT calculates the frequencies in an (maybe) unexpected format
* Data: $[x_1, x_2, ..., x_N]$
* FFT of Data: $[f_0, f_1, ..., f_{N/2}, f_{-N/2}, f_{-N/2+1}, ...,f_{-1}]$
* (for real data negative and positive frequencies are the complex conjugate)

Sometimes we want to have the frequencies centered, which means a `fftshift` operations is needed
* FFT of Data followed by fftshift: $[f_{-N/2}, f_{-N/2+1},..., f_{-1}, f_0, f_1,..., f_{N/2-1}, f_{N/2}]$
"

# ╔═╡ 7276f9e0-80b8-4794-9125-fca43246f818
"""
    my_fft(data::Vector{T})::Vector{Complex{T}}

Calculate the FFT of `data` according to the equation given above.
Note, the output should be a complex Vector!
"""
function my_fft(data::Vector{T})::Vector{Complex{T}} where T
	out = zeros(Complex{T}, size(data))

	# TODO



	
	return out::Vector{Complex{T}}
end

# ╔═╡ e895cff6-4018-495f-9a0e-1e56ddcfd17d
arr = [1f0,2,3,4]

# ╔═╡ 14d68904-4de5-4f1d-81f6-dae8bf50b749
my_fft(arr)

# ╔═╡ 7307dd33-cb79-48c0-9e73-78ed70e9c628
fft(arr)

# ╔═╡ 393d03e2-461d-480a-9719-f33a4eef3a6f
md"### 3.1 Test"

# ╔═╡ df12fd86-f29d-4926-b360-daf4cb40af86
begin
	arr_r = randn((13,))
	PlutoTest.@test fft(arr_r) ≈ my_fft(arr_r)
end

# ╔═╡ 37eedc2e-27bd-4d1c-a83b-6ec14aba55ab
begin
	arr_r2 = randn((14,))
	PlutoTest.@test fft(arr_r2) ≈ my_fft(arr_r2)
end

# ╔═╡ 5313dbed-bfd7-4bb8-b39c-89c55e50d1bb
md"## 3.2 Task
Calculate the FFT of the following dataset, such that center frequency is in the center!

Use `fft` and `fftshift`.
"

# ╔═╡ e1eac3cd-e6ce-4743-9119-5ee0c69b8e00
begin
	data_odd = [1.7142857142857142, 1.643544850921574, 1.0471849694228637, 1.383330948397344, 1.8300917955219322, 1.2951185346328968, 1.086443186817675]
	data_even = [-2.5, 1.5, -2.5, 3.5]
end

# ╔═╡ 481961ac-ad7c-4768-a7aa-eed34143226b
ffts(x) = similar(x)# TODO

# ╔═╡ fc3bb52b-f9a2-42c0-b34b-347d0a137e2d
ffts(data_even)

# ╔═╡ a55851b6-3f8d-4af3-8f09-a3777c56df1b
ffts(data_odd)

# ╔═╡ e46d67be-d974-4f07-8b3e-02fff3e749d0
md"## 3.3 Task

Now do the inverse operation such that `data ≈ iffts(ffts(x))`

Check out `ifftshift`.
"

# ╔═╡ deb9eaf0-ee89-49c5-90b8-570d1a88f8fc
iffts(x) = similar(x) # TODO

# ╔═╡ 2875999a-a61d-491a-8a42-0552a87d3c6a
md"""
### 3.2 and 3.3 Test
"""

# ╔═╡ 333a253e-ce8b-477a-9db2-f99c034ae048
PlutoTest.@test real(ffts(iffts(data_even))) ≈ data_even

# ╔═╡ e81990c8-2194-472c-a102-b0d010f3a266
 PlutoTest.@test real(ffts(iffts(data_odd))) ≈ data_odd

# ╔═╡ f3fc2ea5-033c-4789-8ce6-0c6b15b607b7
 PlutoTest.@test real(iffts(ffts(data_odd))) ≈ data_odd

# ╔═╡ 90596b8a-53d4-457f-9e93-522e6c5dc35d
 PlutoTest.@test real(iffts(ffts(data_even))) ≈ data_even

# ╔═╡ 4e54ddc6-b9cd-472f-b496-18055ee5b667
md"## 3.4 Task
Try to reproduce this $\cos$ curve. 
The reason it looks so _sharp_ is that we only sample slightly above the Nyquist limit.

However, you are only allowed to use `fft`, `fftshift`, `ifftshift`, `ifft` and array manipulation (like `x[34] = 1723.12`).
And you can also take the `real` part at the end.

Please explain with some #comments what you did. 
"

# ╔═╡ 9844ec03-b169-40a0-ad56-1c07b11bb886
begin
	x = range(0, 4, length=21)[begin:end-1]
	y = cos.(x .* 2π)
	plot(x, y, mark="-*")
end

# ╔═╡ 49c39cea-895e-429a-909c-d1438d0dfd00
function reproduce_cos()
	# still wrong
	guess = zeros(20)

	guess[3] = 20
	return real(ift(guess))
end

# ╔═╡ 5c4f40df-cde3-4140-ace2-6c0e1b60ebcf
begin
	plot(x, y, mark="-*")
	plot!(x, reproduce_cos(), mark="-*")
end

# ╔═╡ 79578e4d-7527-459b-8915-ee5cb9ccad24
md"### 3.4 Test"

# ╔═╡ a74c4034-fa19-4360-9f5e-191d5ba4001a
PlutoTest.@test real(reproduce_cos()) ≈ y

# ╔═╡ 8ce249a3-0881-49e5-8ffd-5df21bda50f4
md"## 3.5 Fourier Shift Theorem

The discrete Fourier transform is given by

$$X_k = \sum_{n=0}^{N-1} x_n \exp\left(-i2 \pi \frac{kn}{N} \right)$$
and the inverse by

$$x_n = \frac1{N}\sum_{k=0}^{N-1} X_k \exp\left(i2 \pi \frac{kn}{N} \right)$$.
"

# ╔═╡ 13131381-11b3-4db4-9351-fb7264c94ea4
md"## 3.4 Task
Proof the Fourier shift theorem (with LaTeX) which states

$x_{n+\Delta n} = \mathcal{F}^{-1}\left[ X_{k} \cdot f(k, \Delta n)\right]$

Derive how the function $f(k, \Delta n)$ looks like!

hint: by clicking on the eye symbol next to this task, you can see how to embed LaTeX code in markdown comments.
"

# ╔═╡ bc3ae07f-3877-480f-bd0b-633ce21172e3
md"### Derivation


# TODO
"

# ╔═╡ 8bf6ba5b-1d06-41de-b22a-b58c1eaf2c83
"""
	shift(x::AbstractArray{T, N}, Δn::Number) 

Shifts a (real) array via the Fourier shift theorem by the amount of
`Δn`.

If `Δn ∈ ℕ` then `circshift(x, -Δn) ≈ shift(x, Δn).`
Otherwise a sub-pixel shift is obtained which is equivalent
to a sinc interpolation.


**Hint**: Use `fftfreq` to obtain the frequencies in Fourier space. Otherwise you might run into issues.

## Examples
```julia
julia> shift([1.0, 2.0, 3.0, 4.0], 1)
4-element Vector{Float64}:
 2.0
 3.0
 4.0
 1.0

julia> shift([1.0, 2.0, 3.0, 4.0], -1)
4-element Vector{Float64}:
 4.0
 1.0
 2.0
 3.0

julia> shift([1.0, 2.0, 3.0, 4.0], -0.5)
4-element Vector{Float64}:
 2.5
 1.085786437626905
 2.5
 3.914213562373095
```
"""
function shift(x::AbstractArray{T, N}, Δn::Number) where {T, N}
	# TODO



	
	return similar(x)
end

# ╔═╡ a21d8b13-20db-4df1-a2c7-451f8942466b
begin
	x2 = 0:0.1:3
	y2 = sin.(x2.*3) .+ cos.(x2 .* 4).^2
end;

# ╔═╡ 420d1b15-8dc9-4c9a-b91f-9cc4b8101b92
md"
 $\Delta n$=$(@bind Δn Slider(0:0.01:100, show_value=true))"

# ╔═╡ 54ebf4c4-6692-4a40-a83e-c548440b3b85
begin
	plot(shift(y2, Δn), label="FFT based shift", mark="-*")
	plot!(circshift(y2, round(Int, .-Δn)), label="rounded to integer with circshift", mark="-*")
end

# ╔═╡ bbe52a98-9811-4c3a-9b4e-b90a1f5c163f
md"### 3.5 Test"

# ╔═╡ 73fd7360-ffdf-4870-a301-e39a982e0517
arr4 = randn((37,));

# ╔═╡ 37a815a0-c6cd-4ee8-9f46-ee4aa5259690
PlutoTest.@test circshift(arr4, -12) ≈ shift(arr4, 12)

# ╔═╡ 3c09615a-e9c4-462b-9f23-5d5676aa29b7
PlutoTest.@test FourierTools.shift(arr4, -13.1) ≈ shift(arr4, 13.1)

# ╔═╡ d3f2648a-a736-4812-b721-70298151a1a7
md"# 4 Noise Beyond Bandlimit

As seen in the lecture, noise can move information even beyond
the bandlimit of the optical system.
In this task, we want to demonstrate this behaviour with a simulation.
For that we need a few helper functions, try to fill the missing lines.


## 4.1 rr2 Function 
In the first step, we create a rr2 function.
"

# ╔═╡ 9ef7592b-92f9-48bb-aa3e-9710bf251379
begin
	"""
		rr2(T, s, center=s .÷ 2 + 1)
	
	Returns a 2D array which stores the squared distance to the center.

	## Examples
	```julia
	julia> rr2(Float32, (5,5))
	5×5 Matrix{Float32}:
	 8.0  5.0  4.0  5.0  8.0
	 5.0  2.0  1.0  2.0  5.0
	 4.0  1.0  0.0  1.0  4.0
	 5.0  2.0  1.0  2.0  5.0
	 8.0  5.0  4.0  5.0  8.0
	
	julia> rr2((4,4))
	4×4 Matrix{Float64}:
	 8.0  5.0  4.0  5.0
	 5.0  2.0  1.0  2.0
	 4.0  1.0  0.0  1.0
	 5.0  2.0  1.0  2.0
	
	julia> rr2((3,3), (3,3))
	3×3 Matrix{Float64}:
	 8.0  5.0  4.0
	 5.0  2.0  1.0
	 4.0  1.0  0.0
	```
	"""
	function rr2(T::DataType, s, center=s .÷ 2 .+ 1)
		# TODO


		
		return zeros(s)
	end
	
	function rr2(s, center=s .÷ 2 .+ 1)
		rr2(Float64, s, center)
	end
end

# ╔═╡ b62519f5-26dd-49ca-afa4-2c72de7ef8a6
rr2((4,4))

# ╔═╡ bd45654d-15e0-4717-962b-d5a302de4a9b
PlutoTest.@test rr2(Float32, (4, 3), (2.5, 3)) == Float32[6.25 3.25 2.25; 4.25 1.25 0.25; 4.25 1.25 0.25; 6.25 3.25 2.25]

# ╔═╡ c6ed3889-0d44-4180-9f81-3095055c6e54
PlutoTest.@test rr2((3, 3), (3, 3)) == [8.0 5.0 4.0; 5.0 2.0 1.0; 4.0 1.0 0.0]

# ╔═╡ 5a16e94e-b01b-4897-9631-41d345229be2
md"## 4.2 circ Function
Having a function which creates the radial distance, we now create a function which returns an array which is 1 inside an circle of radius `r` and 0 outside of it.
"

# ╔═╡ 06a666f9-c99e-41e5-b703-a5e826721d04
begin
	"""
		circ(size, radius)
	
	`size` the size of the resulting array
	and `radius` the radius of the circle
	"""
	circ(size, radius) = ones(size) # TODO
end

# ╔═╡ 592e0ed4-bd97-4bb6-a826-b1e2f5ddbd59
circ((5,5), 2)

# ╔═╡ 205a744b-e8a1-4acc-ba6b-7a7307f9a028
md"## 4.2 Test"

# ╔═╡ daf4e474-0db4-4865-859d-96c31e7947f0
PlutoTest.@test circ((10, 2), 2.76) ≈ Bool[0 0; 0 0; 0 0; 1 1; 1 1; 1 1; 1 1; 1 1; 0 0; 0 0]

# ╔═╡ 7fb8fedf-5684-4874-8f17-d01debaff477
PlutoTest.@test circ((4, 5), 2.1) ≈ Bool[0 0 1 0 0; 0 1 1 1 0; 1 1 1 1 1; 0 1 1 1 0]

# ╔═╡ 9b9fb85b-13d0-43a2-ac5f-713593117b96
md" ## 4.3 Frequency Filter
Now we combine all methods together such that we frequency filter
our image.

Procedure:
* Transform image to Fourier space
* take only the inner part of a circle a `radius` of the Fourier space 
* go to real space and take real part
"

# ╔═╡ 2a48072b-fc5a-4044-9aa9-5661483ffbfc
function frequency_filter(img, radius)
	# TODO

	
	return similar(img)
end

# ╔═╡ 643bf329-cce9-4c61-931d-65455b59357d
md"r=$(@bind r Slider(1:512, show_value=true))"

# ╔═╡ 2cafea2d-b4b6-4ede-b899-720a71034fda
gray_show(frequency_filter(img, r))

# ╔═╡ 7467608a-ce6d-4f29-8cd8-d766395f7129
md"## 4.4 Noise Breaks Resolution Limit?
We demonstrate now, that noise beyond the frequency limit still contains information

### Procedure
* First frequency filter our image such that it definetely bandlimited
* Apply Poissonian and Gaussian noise
* Now do the exact opposite of the above low-pass filter and take only the high frequency content!
"

# ╔═╡ 6d7fd369-1f6c-487d-b399-c8cbbfede916
function simulate(img, noise_function, radius)
	# TODO


	
	return similar(img)
end

# ╔═╡ e8bd9cb9-bd20-4619-93da-20e7829a61b4
# anonymous noise function we are going to pass to `simulate`
begin
	p(x) = poisson(x, 100000)
	g(x) = add_gauss(x, 0.1)
end

# ╔═╡ 8a479712-7290-4c00-af7c-44e11ec81d0d
md"r=$(@bind r2 Slider(1:300, show_value=true))"

# ╔═╡ 372e190b-0344-4001-9c72-171610b41e9c
img_noisy = abs.(simulate(img, p, r2));

# ╔═╡ b17cd903-bce6-4291-9d71-ea8a8422295a
gray_show(img_noisy, set_one=true)

# ╔═╡ fc91074a-99fc-4a13-802c-740872bf7502
md"
Mean is $(round(mean(img_noisy), sigdigits=3))

Maximum is $(round(maximum(img_noisy), sigdigits=3))

Minimum is $(round(minimum(img_noisy), sigdigits=3))
"

# ╔═╡ fac6c740-93c1-4e12-9610-f3c621f84043
md"## 5.1 Undersampling
Demonstration of aliasing.

Write a function `undersample` which only takes every `n-th` point of an array, keeping the first pixel of the array.
"

# ╔═╡ 15aa846c-7b81-4094-8b5b-dd6575e82afe
img3 = Float32.(testimage("resolution_test_512"));

# ╔═╡ be3b0811-c05a-4e78-9d34-8c52195dd7d0
"""
	undersample(x, factor)

```julia-repl
julia> undersample([1,2,3,4,5,6,7,8,9,10], 2)
5-element Vector{Int64}:
 1
 3
 5
 7
 9

julia> undersample([1,2,3,4,5,6,7,8,9,10], 3)
4-element Vector{Int64}:
  1
  4
  7
 10

julia> undersample([1,2,3,4,5,6,7,8,9,10], 4)
3-element Vector{Int64}:
 1
 5
 9
```
"""
function undersample(x, factor)
	return ones(size(x) .÷ factor)
end

# ╔═╡ e5188537-b955-4637-bab5-a76b6d74c189
PlutoTest.@test undersample([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 2) == [1, 3, 5, 7, 9]

# ╔═╡ c2b3cdc7-60b5-4287-8268-f95184d452b7
PlutoTest.@test undersample([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 3) == [1, 4, 7, 10]

# ╔═╡ 4e3a272e-0280-4e27-b014-3f9da5ac8039
PlutoTest.@test undersample([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 4) == [1, 5, 9]

# ╔═╡ a76e4081-7c52-4548-b513-ba03954378db


# ╔═╡ f5fbd0da-fe4c-4162-be5e-f9d9d5b27d87
md"## 5.1 Undersampling Visualization

In this plot you can see a nicely oversampled $\sin$ curve.
"

# ╔═╡ a07cb23e-6fca-4c42-a6b6-a525b2b5e232
xs2 = range(0, 10π, length=500)

# ╔═╡ c373046e-04b3-4159-80c0-8fbec141833f
f2(x) = sin(5 * x)

# ╔═╡ 9c1ed9aa-9d52-425d-892d-af3b8a981663
ys2 = f2.(xs2);

# ╔═╡ 75006409-0338-470d-a41d-fa1aba5c4fca
plot(xs2, ys2, xlabel="x pos", ylabel="amplitude")

# ╔═╡ f2b3d8a0-4bea-4443-9bc3-37258106c3d3
md"
Now we only take each `undersample_factor` data point of the series.

`undersample_factor` = $(@bind undersample_factor Slider(1:40, show_value=true))


Analyze what happens once you're below the Nyquist frequency.
See further what happens if you undersample then even more.
"

# ╔═╡ f6b06d6e-b268-45d6-a7e4-d8614cd0884d
y2_undersampled = undersample(ys2, undersample_factor);

# ╔═╡ eab1fe28-f336-4b33-89a3-4c65ed2a818f
plot(fftfreq(length(y2_undersampled))[1:length(y2_undersampled)÷2],abs.(rfft(y2_undersampled))[1:end-1], xlabel="Positive Frequencies", ylabel="Absolute Value of Frequency amplitude", mark="-*")

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
FourierTools = "b18b359b-aebc-45ac-a139-9c0ccbb2871e"
ImageShow = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
Noise = "81d43f40-5267-43b7-ae1c-8b967f377efa"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
TestImages = "5e47fb64-e119-507b-a336-dd2b206d9990"

[compat]
Colors = "~0.12.8"
FFTW = "~1.5.0"
FourierTools = "~0.3.7"
ImageShow = "~0.3.6"
Noise = "~0.3.2"
Plots = "~1.36.2"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.48"
TestImages = "~1.7.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "948be8f4bfb81fb5b88dd00248ba84f862ce6d3f"

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
git-tree-sha1 = "629c6e4a7be8f427d268cebef2a5e3de6c50d462"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.6"

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
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

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

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "aaabba4ce1b7f8a9b34c015053d3b1edf60fa49c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.4.0"

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
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "00a9d4abadc05b9476e937a5557fcce476b9e547"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.69.5"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "bc9f7725571ddb4ab2c4bc74fa397c1c5ad08943"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.69.1+0"

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
git-tree-sha1 = "3f91cd3f56ea48d4d2a75c2a65455c5fc74fa347"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.3"

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
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

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
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

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
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

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
git-tree-sha1 = "5955a002262c08ab155ed2a6f1bb0787fedf5939"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.36.2"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "efc140104e6d0ae3e7e30d56c98c4a927154d684"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.48"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "45f9da1ceee5078267eb273d065e8aa2f2515790"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.3"

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
git-tree-sha1 = "d12e612bba40d189cead6ff857ddb67bd2e6a387"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.1"

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
git-tree-sha1 = "77fea79baa5b22aeda896a8d9c6445a74500a2c2"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.74"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "e59ecc5a41b000fa94423a578d29290c7266fc10"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

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
# ╠═77e6389e-3491-45b6-9f76-459b15a8922d
# ╠═31c31237-c7bd-4fdb-8f7f-219d323f80e7
# ╟─7148e705-8b90-46ea-8763-778897daf9a7
# ╟─83d9c252-2a8d-47a3-b247-a62b6a175201
# ╠═39c623c7-6063-4360-a2b7-95ea42c5262d
# ╟─eacd8fdf-971c-46ca-914f-f73f65e28ca8
# ╠═25dde5cb-9ed6-498e-b511-50ed0640d0c0
# ╠═a6c024a4-7cb2-4f8d-be24-8727d75679da
# ╠═d5c127ee-1ecf-4ed0-9cd9-6fcdefa95847
# ╟─6e95cbcb-0d89-4eca-8599-b266076708a7
# ╠═910cb7fa-b1fc-4106-acc4-3bf56edfb83d
# ╠═cbebb311-c521-4fa6-a909-1a4df8da4d2d
# ╠═8bb737cb-2d59-49ba-9b89-afef6ca7bbf8
# ╠═52bc8b3d-542e-4a35-b38c-35840f4075d4
# ╠═e3448ffe-5fbf-4bed-bf03-6eb9031e6c12
# ╠═d9308666-efaa-4fcc-a083-f3d12dbab1e2
# ╠═493b6bb4-1dfc-4a82-9360-40944b29826c
# ╟─c80a2ce8-6620-4de1-bb78-5af97caf466d
# ╠═e2e470e4-ff73-4dd5-973c-14e1535b01ea
# ╠═0a957f4a-a975-4087-b458-3236940eda43
# ╟─af9866f4-2b03-4dc4-9f8a-e02ae46e2f9c
# ╠═ee4b6b2d-119e-4ac3-a1f5-500ff045bfc5
# ╠═7a481af3-4ec8-4a58-bb84-e9d3d5381570
# ╠═4006ccf7-df8e-440c-8b3c-af525bb50ed4
# ╠═c871cdbd-b560-46f2-918c-45d24d2f9b75
# ╠═f6b8cd92-c1ab-445b-8b4f-04f521b72cad
# ╠═51b2d7be-fe94-44e9-8e67-73aafb636ef1
# ╟─a19c8850-9a7d-4e26-908b-a067a9006e84
# ╟─b0c67bd6-b5ef-4359-8727-1a3bfefcf759
# ╠═7276f9e0-80b8-4794-9125-fca43246f818
# ╠═e895cff6-4018-495f-9a0e-1e56ddcfd17d
# ╠═14d68904-4de5-4f1d-81f6-dae8bf50b749
# ╠═7307dd33-cb79-48c0-9e73-78ed70e9c628
# ╟─393d03e2-461d-480a-9719-f33a4eef3a6f
# ╠═df12fd86-f29d-4926-b360-daf4cb40af86
# ╠═37eedc2e-27bd-4d1c-a83b-6ec14aba55ab
# ╟─5313dbed-bfd7-4bb8-b39c-89c55e50d1bb
# ╠═e1eac3cd-e6ce-4743-9119-5ee0c69b8e00
# ╠═481961ac-ad7c-4768-a7aa-eed34143226b
# ╠═fc3bb52b-f9a2-42c0-b34b-347d0a137e2d
# ╠═a55851b6-3f8d-4af3-8f09-a3777c56df1b
# ╟─e46d67be-d974-4f07-8b3e-02fff3e749d0
# ╠═deb9eaf0-ee89-49c5-90b8-570d1a88f8fc
# ╟─2875999a-a61d-491a-8a42-0552a87d3c6a
# ╠═333a253e-ce8b-477a-9db2-f99c034ae048
# ╠═e81990c8-2194-472c-a102-b0d010f3a266
# ╠═f3fc2ea5-033c-4789-8ce6-0c6b15b607b7
# ╠═90596b8a-53d4-457f-9e93-522e6c5dc35d
# ╟─4e54ddc6-b9cd-472f-b496-18055ee5b667
# ╠═9844ec03-b169-40a0-ad56-1c07b11bb886
# ╠═49c39cea-895e-429a-909c-d1438d0dfd00
# ╠═5c4f40df-cde3-4140-ace2-6c0e1b60ebcf
# ╟─79578e4d-7527-459b-8915-ee5cb9ccad24
# ╠═a74c4034-fa19-4360-9f5e-191d5ba4001a
# ╟─8ce249a3-0881-49e5-8ffd-5df21bda50f4
# ╠═13131381-11b3-4db4-9351-fb7264c94ea4
# ╠═bc3ae07f-3877-480f-bd0b-633ce21172e3
# ╠═8bf6ba5b-1d06-41de-b22a-b58c1eaf2c83
# ╠═a21d8b13-20db-4df1-a2c7-451f8942466b
# ╟─420d1b15-8dc9-4c9a-b91f-9cc4b8101b92
# ╠═54ebf4c4-6692-4a40-a83e-c548440b3b85
# ╟─bbe52a98-9811-4c3a-9b4e-b90a1f5c163f
# ╠═73fd7360-ffdf-4870-a301-e39a982e0517
# ╠═37a815a0-c6cd-4ee8-9f46-ee4aa5259690
# ╠═3c09615a-e9c4-462b-9f23-5d5676aa29b7
# ╟─d3f2648a-a736-4812-b721-70298151a1a7
# ╠═9ef7592b-92f9-48bb-aa3e-9710bf251379
# ╠═b62519f5-26dd-49ca-afa4-2c72de7ef8a6
# ╠═bd45654d-15e0-4717-962b-d5a302de4a9b
# ╠═c6ed3889-0d44-4180-9f81-3095055c6e54
# ╟─5a16e94e-b01b-4897-9631-41d345229be2
# ╠═06a666f9-c99e-41e5-b703-a5e826721d04
# ╠═592e0ed4-bd97-4bb6-a826-b1e2f5ddbd59
# ╟─205a744b-e8a1-4acc-ba6b-7a7307f9a028
# ╠═daf4e474-0db4-4865-859d-96c31e7947f0
# ╠═7fb8fedf-5684-4874-8f17-d01debaff477
# ╟─9b9fb85b-13d0-43a2-ac5f-713593117b96
# ╠═2a48072b-fc5a-4044-9aa9-5661483ffbfc
# ╠═2cafea2d-b4b6-4ede-b899-720a71034fda
# ╠═643bf329-cce9-4c61-931d-65455b59357d
# ╟─7467608a-ce6d-4f29-8cd8-d766395f7129
# ╠═c7150c2c-afeb-4158-b381-3341bb678fd5
# ╠═6d7fd369-1f6c-487d-b399-c8cbbfede916
# ╠═e8bd9cb9-bd20-4619-93da-20e7829a61b4
# ╠═8a479712-7290-4c00-af7c-44e11ec81d0d
# ╠═372e190b-0344-4001-9c72-171610b41e9c
# ╠═b17cd903-bce6-4291-9d71-ea8a8422295a
# ╟─fc91074a-99fc-4a13-802c-740872bf7502
# ╟─fac6c740-93c1-4e12-9610-f3c621f84043
# ╠═15aa846c-7b81-4094-8b5b-dd6575e82afe
# ╠═be3b0811-c05a-4e78-9d34-8c52195dd7d0
# ╠═e5188537-b955-4637-bab5-a76b6d74c189
# ╠═c2b3cdc7-60b5-4287-8268-f95184d452b7
# ╠═4e3a272e-0280-4e27-b014-3f9da5ac8039
# ╠═a76e4081-7c52-4548-b513-ba03954378db
# ╠═f5fbd0da-fe4c-4162-be5e-f9d9d5b27d87
# ╠═a07cb23e-6fca-4c42-a6b6-a525b2b5e232
# ╠═c373046e-04b3-4159-80c0-8fbec141833f
# ╠═9c1ed9aa-9d52-425d-892d-af3b8a981663
# ╠═75006409-0338-470d-a41d-fa1aba5c4fca
# ╠═f2b3d8a0-4bea-4443-9bc3-37258106c3d3
# ╠═f6b06d6e-b268-45d6-a7e4-d8614cd0884d
# ╠═eab1fe28-f336-4b33-89a3-4c65ed2a818f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
