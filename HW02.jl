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

# ╔═╡ a0b4af53-ded5-4d30-9d2b-7c5d1c2bff80
begin
	using ColorSchemes 	# for different colors maps
	using Plots 		# for heatmap
	using TestImages 	# for testimage
	using ImageShow 	# for better rendering of inline images in Pluto
	using Colors # for Gray
	using PlutoUI
end

# ╔═╡ a7007d37-6f30-462b-8be3-70999c980c60
using PlutoTest

# ╔═╡ d9d77790-30d1-4ca2-83ca-d974b10a9431
begin
	using PoissonRandom # for poisson noise (`pois_rand`)
	using Statistics # for mean, var
end

# ╔═╡ 7a0e442c-9f32-4210-9483-63e0dc9fbe49
using FourierTools

# ╔═╡ dd4f1881-4570-4972-bd40-a4123ad4b6b1
# for fitting
using LsqFit

# ╔═╡ 85a15a0e-37d4-4b81-986c-d8a2fe4185c7
md"# Homework 02
In this homework we work on images, colormaps and sensor calibration.
"

# ╔═╡ 9c2160a5-6f62-4ff8-8ecd-88a94226a4fb
"""
    my_show(arr::AbstractArray{<:Real}; set_one=false, set_zero=false)
Displays a real valued array . Brightness encodes magnitude.
Works within Jupyter and Pluto.
## Keyword args
* `set_one=false` divides by the maximum to set maximum to 1
* `set_zero=false` subtracts the minimum to set minimum to 1
"""
function my_show(arr::AbstractArray{<:Real}; set_one=true, set_zero=false)
    arr = set_zero ? arr .- minimum(arr) : arr
    arr = set_one ? arr ./ maximum(arr) : arr
    Gray.(arr)
end

# ╔═╡ 375801bb-7b43-4983-a0ca-211d106ab11c
begin
	# load a image and shift it to negative values as well
	img = Float64.(testimage("moonsurface")) .- 0.5;
	histogram(img[:], xlabel="intensity", ylabel="occurence", title="Histogram Plot", legend=nothing)
end

# ╔═╡ e2768bdd-3346-4f13-88a8-ab8a43094761
md"## 1. Colormaps
Very often in microscopy, we look at grayscale images.
We already learnt that we can use `Gray.(img)` for that.
However, that is only a very simple viewer.
In general, there are much more sophisticated ones like `View5D.jl` or `Napari.jl`.

As you can see below, `Gray` fails pretty quickly and also does not make use of our human color vision.
"

# ╔═╡ cea22309-272b-4696-850a-f49ad3120345
Gray.(img)

# ╔═╡ c66329a6-84f4-4e0e-87a1-8a62f97ba9f7
md" ### 1.1 Heatmap
`img` contains also negative values which `Gray` does not support by default.
A nicer tool is heatmap.
"

# ╔═╡ 450d3dea-cb97-482b-aaab-bb2919e4ee9e
heatmap(img)

# ╔═╡ 41f64a18-bd0d-49a4-8347-36818b45a576
md"## 1.2 Task - Equal Intensity Color Map
We now want to create a RGB color map which always has the same intensity.

Hence, `green_value + red_value + blue_value = const` for all different inputs.

"

# ╔═╡ f6379b6e-7030-4967-afea-13494ccad294
"""
	to_equal_intensity_tuple(value::T) where T

Returns for a given scalar input a tuple of equal intensity.
There is no unique solution, we only require that the sum is always 1.
For example, below would be a valid solution.

## Examples
```julia-repl
julia> to_equal_intensity_tuple(1.0)
(0.04465819873852073, 0.33333333333333326, 0.6220084679281465)

julia> to_equal_intensity_tuple(0.0)
(0.044658198738520505, 0.3333333333333333, 0.6220084679281462)

julia> to_equal_intensity_tuple(0.3f0)
(0.26402956f0, 0.6503522f0, 0.08561832f0)
```
"""
function to_equal_intensity_tuple(value::T) where T
	# with the function signature we can access T (which is the type of `value`)
	# 0 has type Int
	# 1.0 has type Float64
	# 1.0f0 has type Float32
	if value < 0 || value > 1
		error("value must be within [0,1]")
	end

	# fix those three functions such that the property is achieved
	# maybe think about periodic functions ;)
	red_map(x) = x # TODO
	green_map(x) = x # TODO
	blue_map(x) = x # TODO

	return (0,0,1) # TODO, return the correct tuple!
end

# ╔═╡ 2b3c2ab1-9bec-43f5-a008-795ac354e080
md"#### Maybe a visualization helps
"

# ╔═╡ df017bc1-a29c-4bbc-80c3-bc2466cd95f5
begin
	unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))
	x = 0:0.01:1
	r,g,b = unzip(to_equal_intensity_tuple.(x))
	plot(x, r, color=:red, label="red")
	plot!(x, g, color=:green, label="green")
	plot!(x, b, color=:blue, label="blue")
end

# ╔═╡ f05ca6de-b46c-4400-a141-c662a1a16d7c
md"## 1.2 Test"

# ╔═╡ f1d7bf0a-6806-45fb-8213-313416017c18
# cyclic
PlutoTest.@test all(to_equal_intensity_tuple(0f0) .≈ to_equal_intensity_tuple(1f0))

# ╔═╡ 18bb1f0f-5b80-43d5-bc8a-526947a70f2c
PlutoTest.@test to_equal_intensity_tuple(0.1) != to_equal_intensity_tuple(0.2)

# ╔═╡ 76e2f439-1dff-4f2a-9c05-14e10a13fda6
# intensity is 1
PlutoTest.@test 1 ≈ sum(to_equal_intensity_tuple(rand(0f0:1f0)))

# ╔═╡ 5315fbfd-8ade-46d5-9851-f3ee98c349f7
# test for many values
PlutoTest.@test all(1 .≈ sum.(to_equal_intensity_tuple.(rand(0.0:1.0, 100))))

# ╔═╡ 9d4337c8-6cd1-4112-9477-4a40d13fe7e3
md"##### Check if the output is always 0<=output<=1

_Side Note: it took quite a while to engineer this line :D_
"

# ╔═╡ 4ad64115-5b07-4760-963b-ca631970132e
bbb = foldl(&, foldl.(Ref((acc, x) -> acc || 0<=x<=1), to_equal_intensity_tuple.(rand(0.0:1.0, 100)), init=false))

# ╔═╡ c0debce4-d81b-4d5d-8907-53d567cf664a
PlutoTest.@test bbb

# ╔═╡ 4456a207-7f56-4900-95dc-77fe7de31c5d
md"#### Type Stability

The test might be wrong. 
Check what happens with expresions like `2π / 3`. Which type do they have?
How can we convert them to different types?
"

# ╔═╡ a6fbfa34-33c6-4f09-9dcd-773d92a4db0d
PlutoTest.@test NTuple{3, Float32} == typeof(to_equal_intensity_tuple(rand(Float32.(0.0:1.0))))

# ╔═╡ b70be939-7f1d-4e80-9d7f-cc2947d57280
PlutoTest.@test NTuple{3, Float64} == typeof(to_equal_intensity_tuple(rand(Float64.(0.0:1.0))))

# ╔═╡ e9cc6238-5ad9-4a50-93a7-f7065648574a
md"## 1.3 Register Color Map
We are able to create valid outputs, but we still need to pass that to the colormap mechanism.
We need to register the colormap so that we can load it with `heatmap`.

You don't need to know the details, but if you'd like understand, you find those in [ColorSchemes.jl](https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#Make-your-own-colorscheme).
"

# ╔═╡ 633d295e-cc1f-4864-a36b-11905bc47a50
function create_equal_intensity_colormap()
	values = range(0, 1, length=128)
	return map(x -> RGB(x...), to_equal_intensity_tuple.(values))
end

# ╔═╡ 8f0b0ff8-8b9b-4e3a-9b0b-7e10599c5e41
begin
	my_equal_intensity_colormap = create_equal_intensity_colormap()
	loadcolorscheme(:my_equal_intensity_colormap, my_equal_intensity_colormap)
end

# ╔═╡ 6b65dc02-9a84-43c2-9bb5-32104d23f580
heatmap(img, c=:my_equal_intensity_colormap)

# ╔═╡ 160051e5-56b9-413a-a78d-b371a85ff18c
md" ## 1.4 Task
Now we want to create a colormap, which maps values between 0 and 0.5
to blue values and all values above 0.5  to red.
Don't let you confuse by negativ and positive intensities.
`loadcolorschemes` needs the values within the interval [0, 1] and will automatically scale it to the image.

"

# ╔═╡ 1483db2c-6374-4c93-ade4-db8894071ce7
"""
	to_negativ_positive_tuple(value::T) where T


## Examples
```
julia> to_negative_positive_tuple(0.0)
(0.0, 0.0, 1.0)

julia> to_negative_positive_tuple(0.25)
(0.5, 0.5, 1.0)

julia> to_negative_positive_tuple(0.5)
(1.0, 1.0, 1.0)

julia> to_negative_positive_tuple(0.75)
(1.0, 0.5, 0.5)

julia> to_negative_positive_tuple(1.0)
(1.0, 0.0, 0.0)
```
"""
function to_negative_positive_tuple(value::T) where T
	# TODO build the function in analogy to to_equal_intensity_tuple

	# TODO
	# TODO
	
	# TODO
	# TODO
	
	# TODO
	# TODO
	
	# TODO
	# TODO

	return (0, 0, 0) # TODO
end

# ╔═╡ 0afb1220-1fdf-447d-8e79-b1ad7c274dd3
md"
Try to show again a heatmap image. The procedure is exactly the same above, except that you need to exchange the core part of the color generation."

# ╔═╡ 0ff23998-8f11-4a69-a74a-39d3ed73aec8


# ╔═╡ f6e5d3c8-172e-4fb6-89ce-a1827a6ddd6c
begin
	# try to get :my_negative_positive_colormap working
	# see above what to do!

	# TODO
	# TODO
	# TODO
end

# ╔═╡ bd71f8c1-e1e8-4f72-9ad8-8714e9e215f0
heatmap(img, c=:my_negative_positive_colormap)

# ╔═╡ 45d098a5-86e1-45ed-9e27-d59dbc7cfc42
md"## 1.4 Test"

# ╔═╡ a389a422-3946-49ba-ac2a-44ba0ebbb7fe
PlutoTest.@test to_negative_positive_tuple(0.0) == (0.0, 0.0, 1)

# ╔═╡ cf832e74-1029-4554-94ff-876e7ce622ab
PlutoTest.@test to_negative_positive_tuple(1.0) == (1.0, 0.0, 0.0)

# ╔═╡ 968529e7-c30c-4ba9-96cc-c3468bbbf7e3
PlutoTest.@test to_negative_positive_tuple(0.5) == (1.0, 1.0, 1.0)

# ╔═╡ 65aab672-9fc2-48f9-b669-e6f5299b1e99
PlutoTest.@test to_negative_positive_tuple(0.75) == (1.0, 0.5, 0.5)

# ╔═╡ 00ee1b88-41fc-41d6-b6a5-91a68560384d
PlutoTest.@test to_negative_positive_tuple(0.25) == (0.5, 0.5, 1)

# ╔═╡ b7e01e25-4b4c-4021-ac52-b5fa4d54c1fa
PlutoTest.@test typeof(to_negative_positive_tuple(0.25f0)) == Tuple{Float32, Float32, Float32}

# ╔═╡ ac0c999b-afb3-4b0c-b981-1bde3275e182
PlutoTest.@test typeof(to_negative_positive_tuple(0.0)) == Tuple{Float64, Float64, Float64}

# ╔═╡ 279cd4f9-9d0e-4bef-80c9-907652371f12
md"## 2 Sensor Simulation

In this part, we want to simulate a sensor. Later, we try to calibrate our artifical sensor.

In short: Each pixel in the sensor measures a number of photons (`n_photons`) which is affected by Poisson noise.
Afterwards, the measured value is altered by additive Gauss noise with a certain
standard deviation `σ`.
This value is then converted to ADU units (analog to digital units) with a certain `gain` (linear conversion factor).
Additionally, each sensor has an `offset` which is finally added.
"

# ╔═╡ ffbc76ee-5c5b-460f-a7dc-ac50500dede9
md"## 2.1 Task
Fill in the missing parts.
"

# ╔═╡ 0b6561db-6286-49b6-bbfb-98f56cbafd1b
"""
      simulate_pixel(n_photons; read_σ=5, offset=0.5, gain=0.1)

In this function we simulate a single pixel.
It transforms a photon number to a digital unit (ADU).

The input is a integer number of photons which is then altered with poisson noise (use `pois_rand` for that).
This number is changed with additive Gauss noise (see HW01).
Finally, the result of this operation is multiplied by the linear `gain`.
At the end, we add `offset`.



* `n_photons`: number of photons
* `read_σ`: additive read noise of the sensor (applied before `gain`)
* `gain`: how much digital output per photon
* `offset`: how much digital offset


## Example Outputs - can vary in your case!
```julia
julia> simulate_pixel(10)
0.8570431305838304

julia> simulate_pixel(10)
2.060753446522642

julia> simulate_pixel(10)
1.020461460534897

julia> simulate_pixel(100)
8.96276245560446

julia> simulate_pixel.(ones((3,3)))
3×3 Matrix{Float64}:
 0.396078  0.96849   -0.0387851
 1.09982   0.210089   0.685605
 0.606973  1.42573   -0.675297
```
"""
function simulate_pixel(n_photons; read_σ=5, offset=0.5, gain=0.1)
	# try to compose the different operations
	# TODO
	# TODO
	return 0 # TODO
end

# ╔═╡ 3b1bbda7-dfa0-4181-b0f7-fc920f166a4b
single_output = simulate_pixel(100) 

# ╔═╡ 1e89c538-c8a4-4229-b026-96bfe15f1eba


# ╔═╡ eac17ab8-9ec4-445e-b317-e2b9dcaea4df


# ╔═╡ 5df0c449-c46e-4719-bff5-d42220991e8f
md"## 2.1 Test"

# ╔═╡ 02b3bafb-fc03-42ff-a281-47216a04f6ce
# check offset
PlutoTest.@test 0.0123 ≈ simulate_pixel(0, read_σ=0, offset=0.0123, gain=1)

# ╔═╡ 1d3e3c6f-7920-4b2c-a432-105fa1f131f7
# check mean poisson
PlutoTest.@test 123 ≈ mean(simulate_pixel.(ones((1000,)) * 123, read_σ=0, offset=0.0123, gain=1)) rtol=0.1

# ╔═╡ b01359d9-c1b2-455f-897c-8ab6bed4f20f
# check variance poisson
PlutoTest.@test 123 ≈ var(simulate_pixel.(ones((1000,)) * 123, read_σ=0, offset=0.0123, gain=1)) rtol=0.1

# ╔═╡ a598763f-63f9-48c9-9e5b-ba6256775221
# check readnoise
PlutoTest.@test 12.3 ≈ std(simulate_pixel.(zeros((1000,)) * 123, read_σ=12.3, offset=0, gain=1)) rtol=0.1

# ╔═╡ 93567e80-181f-4629-8a27-2f37e5e20b62
# check gain
PlutoTest.@test 42 * 12.3 ≈ std(simulate_pixel.(zeros((1000,)) * 123, read_σ=12.3, offset=0, gain=42)) rtol=0.1

# ╔═╡ 9ad71b79-872d-4c96-b4eb-3c2e841f295d
md"## 2.2 Variance Mean Projection
Now we want to do a Variance Mean Projection to fit the values of `offset` and `gain` since we often don't know them for real sensors.

For this test, you take a sample and in practice you defocus your microscope heavily.
Via that trick, you have an image which covers all values between very dark and very bright.
Hence, we expect that our image sensors both measures a large dynamic range.
"

# ╔═╡ 70fa1c62-b5df-49dd-ae66-b49cb9bf6c0a
# that's already done for you
begin
	ideal = resample(Float32.(testimage("resolution_test_512"))[50:450, 50:450], (150, 150)) # random amount of photons in the interval 0:500
	x_cords = range(-10, 10, length=size(ideal, 1))
	gauss = Float32.(exp.(-(x_cords.^2 .+ x_cords'.^2)))
	gauss ./= sum(gauss)
	img_blurry = conv_psf(ideal, gauss)
	img_blurry .-= minimum(img_blurry)
	img_blurry ./= maximum(img_blurry) 
	img_blurry = round.(Ref(Int), img_blurry .* 1000) # round to integers
	Gray.(img_blurry ./ 500)
end;

# ╔═╡ 1e3f86d7-a952-41fc-9f67-c8ce87d14315
md"###### That's an ideal image, where we have both very bright and very dark regions"

# ╔═╡ 120846a5-7892-4f79-8e95-bde9bf06df25
# our blurry sample which emits 1000 photons in maximum
my_show(img_blurry)

# ╔═╡ 3a6757db-050c-4923-a6bd-5c7c33581a93
histogram(img_blurry[:], xlabel="intensity", ylabel="occurence", title="Histogram Plot", legend=nothing)

# ╔═╡ ee7b9e26-c81d-4e09-9553-e9cded4db31d
md"## 2.2 Task - Apply to image
Now apply the function `simulate_pixel` to `img_blurry`.
Use the automatic broadcasting of Julia for that.
"

# ╔═╡ 6fddf2a4-45f0-4243-aa0c-df5ce2c354a9
simulated_img = similar(img) # TODO, fix that line such that it uses `simulate_pixel`. Use broadcasting mechanism

# ╔═╡ 48fb4f01-a3fa-4c7b-abb0-6e6f2b63dac7
md"Apply the same trick (broadcasting) for a 2D image to this 3D image"

# ╔═╡ f9b327ef-e184-4e1f-af3c-184f408f04bd
md"## 2.3 Task - Calculate mean and variance
Calculate the mean and variance of those images, but only along the third dimension!
Either use a library (Statistics) or write it yourself :)
"

# ╔═╡ 087bd202-d950-4573-a010-72f9326bd410
 arr_test = [1; 2;;; 3; 2]

# ╔═╡ 392c1a92-b297-4f20-a5ff-edbfe1f124f1
md"
In this case, the mean along the third dimension would be"

# ╔═╡ ea1687f3-5cc2-4dcb-9317-019c33d57291
[2.0; 2.0;;;]

# ╔═╡ ec0289d5-427c-4e19-b63a-361ed7165d42
md"#### Finally, we produce a scatter plot"

# ╔═╡ 1df9763e-81be-4f24-b439-e01559d78a9b
md" ## 2.4 Fitting
Having our noisy dataset, we want to extract the `gain` of the detector since the gain translate from a photon count to an eletrical signal (ADU).

Why does plotting the variance $\sigma^2$ over the mean $\mu$ allows to fit the gain?
We know, that in a Poisson distribution the expected mean $\mu$ is always equal to the variance $\sigma^2$.

If we now take 10 captures of the same image, we always expect the same mean $\mu$ for a certain pixel. Of course, those are 10 evaluation of a random measurement process, therefore we observe fluctuations. Hence, for each pixel we can calculate the mean and the variance. 

Since the detector multiplies the measured photons with a certain gain, we do not measure the photons but the digital value.
Hence the measured mean is

$\hat \mu = \text{gain} \cdot \mu$

where $\mu$ would be the expected photon number.

However, the measured variance is 

$$\hat \sigma^2 = \sum_{i} (\hat \mu - \text{gain} \cdot x_i)^2 = \text{gain}^2 \sigma^2$$
where $\sigma^2$ would be the variance of the photon number.

The slope of our variance over mean is then 

$$\frac{\Delta \hat \sigma^2}{\Delta \hat \mu} = \text{gain}$$

Dividing $\hat \mu$ by $\text{gain}$ returns the measured photon number.

"

# ╔═╡ d6ea97c0-42a2-4119-94cb-c59d07c84992
md" ## 2.5 LsqFit.jl
We now want to fit a linear function with a slope and a offset on our data.
Try to find out [here](https://github.com/JuliaNLSolvers/LsqFit.jl) how that works
"

# ╔═╡ 529cf3ab-090d-4ed3-b236-eed577554bbb
@. model(x, p) = nothing # TODO

# ╔═╡ ca78fb0c-3cd6-461e-a4c7-5bbd61589001
p0 = nothing # TODO (maybe you need to convert it to Float32, maybe not)

# ╔═╡ 709c3594-d979-4a09-876f-7bbf4ad4aa18
fit = nothing # TODO

# ╔═╡ b2b821b3-e95b-48ec-aee0-acea4e6b01dc
md"##  Final Check
If you slide those parameters, you change the simulated values.
"

# ╔═╡ e68aa618-f5b0-4834-b492-268962b4c5be
md"
`offset`=$(@bind offset Slider(0.5:1:100, show_value=true))

`read_σ`=$(@bind read_σ Slider(2:1:20, show_value=true))

`gain`=$(@bind gain Slider(0.01:0.01:1, show_value=true))
"

# ╔═╡ 1140ab5b-6e3a-4fac-baf0-a150a07d7293
begin
	imgs = repeat(img_blurry, outer=(1, 1, 10));
	images = simulate_pixel.(imgs; offset, read_σ, gain);
end;

# ╔═╡ 1e97f466-796e-4837-b0bf-2a3c4edd89c0
md" In practice we would need to take a few images (like 10)

Hence we repeat the image 10 times digitally
 
The output is an array with `size(imgs) ==` $(size(imgs))
"

# ╔═╡ 3a7737ac-7884-41ba-863c-b352dceb73c9
begin
	# mean
	μ_arr = randn(size(imgs)[1:2]) # calculate the mean of `images` along dim 3
	μ = μ_arr[:]
end

# ╔═╡ 501a5a71-9cff-4338-b575-0f0be06d064f
begin
	# variance
	variance_arr = randn(size(imgs)[1:2])# calculate the variance of `images` along dim 3
	variance = variance_arr[:]
end

# ╔═╡ 41bdbfb9-06b9-4f4b-bc97-31ed2e6ae695
scatter(μ, variance, xlabel="mean", ylabel="variance")

# ╔═╡ 3ffb57c7-825b-4213-a847-7cc4e9466875
begin
	scatter(μ, variance, label="data")
	plot!(μ, model(μ, fit.param), linewidth=5, label="fit")
end

# ╔═╡ 0adf654e-f238-432c-ace4-50ad5d896870
md"The final solution is gain = $(round(fit.param[1], sigdigits=4)) ± $(round(stderror(fit)[1], sigdigits=2))"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ColorSchemes = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
FourierTools = "b18b359b-aebc-45ac-a139-9c0ccbb2871e"
ImageShow = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PoissonRandom = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
TestImages = "5e47fb64-e119-507b-a336-dd2b206d9990"

[compat]
ColorSchemes = "~3.19.0"
Colors = "~0.12.8"
FourierTools = "~0.3.6"
ImageShow = "~0.3.6"
LsqFit = "~0.13.0"
Plots = "~1.35.8"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.48"
PoissonRandom = "~0.4.3"
TestImages = "~1.7.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "aae49aa9ea12e8abe660eb3136fd9ac60c0f38c7"

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
git-tree-sha1 = "732cddf5c7a3d4e7d4829012042221a724a30674"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.24"

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
git-tree-sha1 = "6d8c6e65b3a57302f094c8d431667b6a4cd2e4cb"
uuid = "26cce99e-4866-4b6d-ab74-862489e035e0"
version = "0.6.6"

[[deps.BitFlags]]
git-tree-sha1 = "84259bb6172806304b9101094a7cc4bc6f56dbc6"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.5"

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

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

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

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "3ca828fe1b75fa84b021a7860bd039eaea84d2f2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.3.0"

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
git-tree-sha1 = "46d2680e618f8abd007bce0c3026cb0c4a8f2032"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.12.0"

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

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "8b7a4d23e22f5d44883671da70865ca98f2ebf9d"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.12.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "04db820ebcfc1e053bd8cbb8d8bccf0ff3ead3f7"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.76"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "c36550cb29cbe373e95b3f40486b9a4148f89ffd"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.2"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

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
git-tree-sha1 = "5a2cff9b6b77b33b89f3d97a4d367747adce647e"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.15.0"

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
git-tree-sha1 = "187198a4ed8ccd7b5d99c41b69c679269ea2b2d4"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.32"

[[deps.FourierTools]]
deps = ["ChainRulesCore", "FFTW", "IndexFunArrays", "LinearAlgebra", "NDTools", "NFFT", "PaddedViews", "ShiftedArrays"]
git-tree-sha1 = "c4b970babd62198d43fc1a7128937a16142dedb2"
uuid = "b18b359b-aebc-45ac-a139-9c0ccbb2871e"
version = "0.3.6"

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
git-tree-sha1 = "a97d47758e933cd5fe5ea181d178936a9fc60427"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.5.1"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

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
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "OptimBase", "Random", "StatsBase"]
git-tree-sha1 = "00f475f85c50584b12268675072663dfed5594b2"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.13.0"

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
git-tree-sha1 = "43f9be9c281179fe44205e2dc19f22e71e022d41"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.15"

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
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

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
git-tree-sha1 = "5c0b972bf8700cfea7fc41b2ad5ca782c356366f"
uuid = "efe261a4-0d2b-5849-be55-fc731d526b0d"
version = "0.13.2"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

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
git-tree-sha1 = "3c3c4a401d267b04942545b1e964a20279587fd7"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.0"

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

[[deps.OptimBase]]
deps = ["NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "9cb1fee807b599b5f803809e85c81b582d2009d6"
uuid = "87e2bd06-a317-5318-96d9-3ecbac512eee"
version = "2.0.2"

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

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

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
deps = ["Dates"]
git-tree-sha1 = "6c01a9b494f6d2a9fc180a08b182fcb06f0958a0"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.2"

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
git-tree-sha1 = "b434dce10c0290ab22cb941a9d72c470f304c71d"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.35.8"

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
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "97aa253e65b784fd13e83774cadc95b38011d734"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.6.0"

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
git-tree-sha1 = "249df6fb3520492092ccebe921829920215ab205"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.9"

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

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

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
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

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
git-tree-sha1 = "f86b3a049e5d05227b10e15dbb315c5b90f14988"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.9"

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

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5783b877201a82fc0014cbf381e7e6eb130473a4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.0.1"

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
git-tree-sha1 = "70e6d2da9210371c927176cb7a56d41ef1260db7"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.1"

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
# ╟─85a15a0e-37d4-4b81-986c-d8a2fe4185c7
# ╠═a0b4af53-ded5-4d30-9d2b-7c5d1c2bff80
# ╠═9c2160a5-6f62-4ff8-8ecd-88a94226a4fb
# ╠═375801bb-7b43-4983-a0ca-211d106ab11c
# ╟─e2768bdd-3346-4f13-88a8-ab8a43094761
# ╠═cea22309-272b-4696-850a-f49ad3120345
# ╟─c66329a6-84f4-4e0e-87a1-8a62f97ba9f7
# ╠═450d3dea-cb97-482b-aaab-bb2919e4ee9e
# ╟─41f64a18-bd0d-49a4-8347-36818b45a576
# ╠═f6379b6e-7030-4967-afea-13494ccad294
# ╟─2b3c2ab1-9bec-43f5-a008-795ac354e080
# ╠═df017bc1-a29c-4bbc-80c3-bc2466cd95f5
# ╟─f05ca6de-b46c-4400-a141-c662a1a16d7c
# ╠═a7007d37-6f30-462b-8be3-70999c980c60
# ╠═f1d7bf0a-6806-45fb-8213-313416017c18
# ╠═18bb1f0f-5b80-43d5-bc8a-526947a70f2c
# ╠═76e2f439-1dff-4f2a-9c05-14e10a13fda6
# ╠═5315fbfd-8ade-46d5-9851-f3ee98c349f7
# ╟─9d4337c8-6cd1-4112-9477-4a40d13fe7e3
# ╠═4ad64115-5b07-4760-963b-ca631970132e
# ╠═c0debce4-d81b-4d5d-8907-53d567cf664a
# ╟─4456a207-7f56-4900-95dc-77fe7de31c5d
# ╠═a6fbfa34-33c6-4f09-9dcd-773d92a4db0d
# ╠═b70be939-7f1d-4e80-9d7f-cc2947d57280
# ╟─e9cc6238-5ad9-4a50-93a7-f7065648574a
# ╠═633d295e-cc1f-4864-a36b-11905bc47a50
# ╠═8f0b0ff8-8b9b-4e3a-9b0b-7e10599c5e41
# ╠═6b65dc02-9a84-43c2-9bb5-32104d23f580
# ╟─160051e5-56b9-413a-a78d-b371a85ff18c
# ╠═1483db2c-6374-4c93-ade4-db8894071ce7
# ╟─0afb1220-1fdf-447d-8e79-b1ad7c274dd3
# ╠═0ff23998-8f11-4a69-a74a-39d3ed73aec8
# ╠═f6e5d3c8-172e-4fb6-89ce-a1827a6ddd6c
# ╠═bd71f8c1-e1e8-4f72-9ad8-8714e9e215f0
# ╟─45d098a5-86e1-45ed-9e27-d59dbc7cfc42
# ╠═a389a422-3946-49ba-ac2a-44ba0ebbb7fe
# ╠═cf832e74-1029-4554-94ff-876e7ce622ab
# ╠═968529e7-c30c-4ba9-96cc-c3468bbbf7e3
# ╠═65aab672-9fc2-48f9-b669-e6f5299b1e99
# ╠═00ee1b88-41fc-41d6-b6a5-91a68560384d
# ╠═b7e01e25-4b4c-4021-ac52-b5fa4d54c1fa
# ╠═ac0c999b-afb3-4b0c-b981-1bde3275e182
# ╟─279cd4f9-9d0e-4bef-80c9-907652371f12
# ╠═d9d77790-30d1-4ca2-83ca-d974b10a9431
# ╟─ffbc76ee-5c5b-460f-a7dc-ac50500dede9
# ╠═0b6561db-6286-49b6-bbfb-98f56cbafd1b
# ╠═3b1bbda7-dfa0-4181-b0f7-fc920f166a4b
# ╠═1e89c538-c8a4-4229-b026-96bfe15f1eba
# ╠═eac17ab8-9ec4-445e-b317-e2b9dcaea4df
# ╟─5df0c449-c46e-4719-bff5-d42220991e8f
# ╠═02b3bafb-fc03-42ff-a281-47216a04f6ce
# ╠═1d3e3c6f-7920-4b2c-a432-105fa1f131f7
# ╠═b01359d9-c1b2-455f-897c-8ab6bed4f20f
# ╠═a598763f-63f9-48c9-9e5b-ba6256775221
# ╠═93567e80-181f-4629-8a27-2f37e5e20b62
# ╟─9ad71b79-872d-4c96-b4eb-3c2e841f295d
# ╠═7a0e442c-9f32-4210-9483-63e0dc9fbe49
# ╠═70fa1c62-b5df-49dd-ae66-b49cb9bf6c0a
# ╟─1e3f86d7-a952-41fc-9f67-c8ce87d14315
# ╠═120846a5-7892-4f79-8e95-bde9bf06df25
# ╠═3a6757db-050c-4923-a6bd-5c7c33581a93
# ╟─ee7b9e26-c81d-4e09-9553-e9cded4db31d
# ╠═6fddf2a4-45f0-4243-aa0c-df5ce2c354a9
# ╟─1e97f466-796e-4837-b0bf-2a3c4edd89c0
# ╟─48fb4f01-a3fa-4c7b-abb0-6e6f2b63dac7
# ╠═1140ab5b-6e3a-4fac-baf0-a150a07d7293
# ╟─f9b327ef-e184-4e1f-af3c-184f408f04bd
# ╠═087bd202-d950-4573-a010-72f9326bd410
# ╟─392c1a92-b297-4f20-a5ff-edbfe1f124f1
# ╠═ea1687f3-5cc2-4dcb-9317-019c33d57291
# ╠═3a7737ac-7884-41ba-863c-b352dceb73c9
# ╠═501a5a71-9cff-4338-b575-0f0be06d064f
# ╟─ec0289d5-427c-4e19-b63a-361ed7165d42
# ╠═41bdbfb9-06b9-4f4b-bc97-31ed2e6ae695
# ╟─1df9763e-81be-4f24-b439-e01559d78a9b
# ╠═dd4f1881-4570-4972-bd40-a4123ad4b6b1
# ╟─d6ea97c0-42a2-4119-94cb-c59d07c84992
# ╠═529cf3ab-090d-4ed3-b236-eed577554bbb
# ╠═ca78fb0c-3cd6-461e-a4c7-5bbd61589001
# ╠═709c3594-d979-4a09-876f-7bbf4ad4aa18
# ╠═3ffb57c7-825b-4213-a847-7cc4e9466875
# ╟─b2b821b3-e95b-48ec-aee0-acea4e6b01dc
# ╟─e68aa618-f5b0-4834-b492-268962b4c5be
# ╠═0adf654e-f238-432c-ace4-50ad5d896870
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
