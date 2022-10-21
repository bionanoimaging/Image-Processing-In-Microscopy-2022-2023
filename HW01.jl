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

# ╔═╡ 8a11c2b7-b6e5-4c26-9e06-ad3ab958db9c
using TestImages, ImageShow, Colors, Statistics, PlutoTest, PlutoUI

# ╔═╡ 3ccbf4fb-0707-4ce8-83d8-b1af555f3fe1
using PoissonRandom

# ╔═╡ 62c6e2f0-0007-40ec-81c5-61cda80f59e5
using ImageFiltering, IndexFunArrays

# ╔═╡ 31cc54e9-d47b-4df7-a785-dad0120e8cdb
md"## Load Packages"

# ╔═╡ 50044e20-1dda-444f-9011-d9a55e182e11


# ╔═╡ b731d80c-3012-49f5-8e63-9b7cea0ed667
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

# ╔═╡ d8c83f2c-76ea-42b9-bd39-7fa5cafcd4d3
md"# Homework 01"

# ╔═╡ 334ebdc5-dc11-4e83-ad12-a1961420a97b
begin
	img = Float64.(testimage("mandril_gray"))
	my_show(img)
end

# ╔═╡ 2122c12f-0832-4f19-9a77-3047d5311cba
md"## 1. Add Noise To Images

This homework is dedicated to introduce into the basics of Julia.
As examples, we want to apply and remove noise from images
"

# ╔═╡ a7d7a415-7464-43fb-a2bf-8e3f2d999c32
md"## Task 1.1  - Gaussian Noise
In the first step, we want to add Gaussian noise with a certain 
standard deviation $\sigma$. Try to get yourself used to the Julia documentation
and check out whether such a function exists already and how you can apply it
to your image/array.

Codewords: **Julia, random, random normal distribution**
"

# ╔═╡ de104b0a-28a4-42de-ae8d-fa857a8f32e0
"""
	add_gauss_noise(img, σ=1)

This function adds normal distributed noise to `img`.
`σ` is an optional argument
"""
function add_gauss_noise(img, σ=one(eltype(img)))
	# todo
	return similar(img)
end

# ╔═╡ f89bc1e4-d406-4550-988b-71496b64035a
md"
Now we want to use a for loop inside this function.
Try to find out how to write for loops to iterate through an array.
"

# ╔═╡ c0b5ff08-7342-4c25-9791-be3a4918c400
"""
	add_gauss_noise_fl!(img, σ=1)

This function adds normal distributed noise to `img`.
`σ` is an optional argument.
This function is memory efficient by using for loops.


`!` means that the input (`img`) is modified.
Therefore, don't return a new array but instead modify the existing one!
The bang (!) is a convention in Julia that a function modifies the input.
"""
function add_gauss_noise_fl!(img, σ=1)
	# todo
	return img
end

# ╔═╡ d56bc500-f5cf-40f3-8c0b-b4ae26a58367
md"### Tests 1.1
Don't modify but rather take those tests as input whether you are correct.
Green is excellent, red not :(
"

# ╔═╡ 1116455e-4698-489a-add2-4e1d8ecbb546
my_show(add_gauss_noise(img, 0.15), set_one=true)

# ╔═╡ cc7bbb55-b8dd-454b-b170-45137ab0c0b5
md"###### Check wether mean and standard deviation are correct"

# ╔═╡ c997c248-fba8-4ed8-ab75-179d248e6b83
PlutoTest.@test ≈(0.3, std(add_gauss_noise(ones((512, 512)), 0.3)), rtol=0.1)

# ╔═╡ 07bc00b7-81ad-4645-8fbe-8545ccdfe7e4
PlutoTest.@test ≈(0.3, std(add_gauss_noise_fl!(ones((512, 512)), 0.3)), rtol=0.1)

# ╔═╡ d317ed1a-1867-4293-8db4-51352212929c
PlutoTest.@test ≈(1, mean(add_gauss_noise(ones((512, 512)), 0.3)), rtol=0.1)

# ╔═╡ e69a8b77-833a-4498-8bc9-e4c970f6529d
md"## Task 1.2 - Poisson Noise
In microscopy and low light imaging situation, the dominant noise term is usually Poisson noise which we want to simulate here.

For adding Poisson noise we use: [PoissonRandom.jl](https://github.com/SciML/PoissonRandom.jl)

Read the documentation of this package how to generate Poisson Random numbers
"

# ╔═╡ 1e4b32cb-00ea-439f-900e-abadc850be95
"""
	add_poisson_noise!(img, scale_to=nothing)

This function adds poisson distributed noise to `img`.

Before adding noise, it scales the maximum value to `scale_to` and 
divides by it afterwards.
With that we can set the number of events (like a photon count).

Differently said: the `scale_to` applies Noise to an array equivalent to the noise level of an array with maximum peak `scale_to`.

If `isnothing(scale_to) == true`, we don't modify/scale the array.

`!` means that the input is modified.
"""
function add_poisson_noise!(img, scale_to=nothing)
	# todo
	return img
end

# ╔═╡ 6cf8ff2a-4255-4c53-9376-5a0d2d372569
add_poisson_noise!(10 .* ones((3, 3)))

# ╔═╡ 9a3feb13-0424-448c-934b-45c476cb0096
add_poisson_noise!(10 .* ones((3, 3)), 100)

# ╔═╡ e3fcd381-9652-4e1d-8585-a91f87b73ce5
md"### Test 1.2 - Poisson Noise

"

# ╔═╡ f3efcfcf-4516-43d2-9375-d8863b4b6f5b
my_show(add_poisson_noise!(100 .* img), set_one=true)

# ╔═╡ 5b68030f-557e-449e-8baf-cf7ae65e4425
mean(add_poisson_noise!(150 .* ones(Float64, (512, 512))))

# ╔═╡ 78a39a77-bc9d-427a-a113-4e64b0c9d311
PlutoTest.@test ≈(150, mean(add_poisson_noise!(150 .* ones(Float64, (512, 512)))), rtol=0.05)

# ╔═╡ 41118b3c-3e49-4fc7-8c80-7906b7cf91fe
PlutoTest.@test ≈(√(150), std(add_poisson_noise!(150 * ones(Float64, (512, 512)))), rtol=0.05)

# ╔═╡ f8be3a4b-d8f0-40cb-869e-81c17e58327a
md"###### Consider those tests as bonus"

# ╔═╡ 65a8981d-890b-4590-a913-6890ecf3817d
PlutoTest.@test ≈(150, 150 * mean(add_poisson_noise!(ones((512, 512)), 150)), rtol=0.05)

# ╔═╡ 272c8afb-667f-49c7-af48-2feb1b9d06f7
PlutoTest.@test ≈(√(150), 150 * std(add_poisson_noise!(ones((512, 512)), 150)), rtol=0.05)

# ╔═╡ 59482e49-d784-4303-9487-bf37f6a0462e
md"## 1.3 Hot Pixels

Another issue are hot pixels which result in a pixel with maximum value. This can be due to damaged pixels or some other noise (radioactivity, ...). Often this is called Salt (because of the maximum brightness) noise.
" 

# ╔═╡ 6e050fc2-4db2-4e6e-abd4-2812b47c070f
"""
	add_hot_pixels!(img, probability=0.1; max_value=one(eltype(img)))

Add randomly hot pixels. The probability for each pixel to be hot,
should be specified by `probability`.
`max_value` is a keyword argument which is the value the _hot_ pixel will have.
"""
function add_hot_pixels!(img, probability=0.1; max_value=1)
	# todo
	return img
end

# ╔═╡ 22985061-9740-45e9-af56-e0787dadba7b
my_show(add_hot_pixels!(copy(img)))

# ╔═╡ b6c25124-3165-4c85-9b1c-da68bd68b59f
my_show(add_hot_pixels!(copy(img), 0.5; max_value=10))

# ╔═╡ 76e83c0b-e9a8-46bd-8a7e-cac39ccaefff
my_show(add_hot_pixels!(copy(img), 0.9))

# ╔═╡ 42f92a8c-54bb-4033-9ef2-f5ffb78f4f3a
md"### 1.3 Test"

# ╔═╡ a9ad9d52-53c9-44ce-8edc-af8d7bfdc81f
PlutoTest.@test sum(map(x -> x ≈ 2, add_hot_pixels!(copy(img), 0.5, max_value=2))) ≈ length(img) .* 0.5 rtol=0.05

# ╔═╡ 373a3e39-fd8c-4cea-88bb-2b9c3734f14f
md"## 2. Remove Noise From Images"

# ╔═╡ 2c274829-ef33-4f34-9ba0-63d3b8cd7343
md"## 2.1 Remove Noise with Gaussian Blur Kernel

One option to remove noise is, is to blur the image with a Gaussian kernel.
We can convolve a small gaussian (odd sized) Kernel over the array.
For that we are using `ImageFiltering.imfilter` for the convolution.

To create a Gaussian shaped function, checkout `IndexFunArrays.normal`. 
You probably want to normalize the sum of it again
"

# ╔═╡ 52ebebfc-940b-47e2-a059-40ec996f2f98
begin
	img_g = add_gauss_noise(img, 0.1)
	img_p = add_poisson_noise!(100 .* copy(img)) ./ 100
	img_h = add_hot_pixels!(copy(img))
end;

# ╔═╡ 2d6f25a5-b902-4e46-ac3f-cf452f525912
my_show([img_g img_p img_h])

# ╔═╡ 75aa8212-8ba3-4fba-b8b2-af2dd9c82adf
# check out the help
imfilter

# ╔═╡ ded0b153-3ba6-4ad6-906b-87930e2b82d0
# check out the help
normal

# ╔═╡ 6acf7ba6-445d-4ac0-a83f-6916318c783c
function gaussian_noise_remove(arr; kernel_size=(3,3), σ=1)
	# todo
	return similar(arr)
end

# ╔═╡ 4f683117-9a09-4837-b0af-5818ac2e62fd
md"σ = $(@bind σ Slider(0.01:0.1:10, show_value=true))"

# ╔═╡ 01203c2c-4e67-492c-9d9c-dc84a84afe15
md"kernel size = $(@bind ks Slider(3:2:21, show_value=true))"

# ╔═╡ daae9a7d-99f4-4d34-a7fd-886fa754d59b
begin
	img_p_gauss = gaussian_noise_remove(img_p, kernel_size=(ks,ks), σ=σ);
	img_g_gauss = gaussian_noise_remove(img_g, kernel_size=(ks,ks), σ=σ);
	img_h_gauss = gaussian_noise_remove(img_h, kernel_size=(ks,ks), σ=σ);
end;

# ╔═╡ 92d43328-29a2-45b2-ae71-e169dec02316
my_show([img_p_gauss img_g_gauss img_h_gauss])

# ╔═╡ 512486eb-1919-446a-9f24-e9713546d237
md"### 2.1 Test"

# ╔═╡ 9b83d92d-f838-4329-b0b1-62bf7a18f123
arr_rand = add_gauss_noise(ones((500, 500)), 0.2);

# ╔═╡ bbbed1b7-a406-44a7-8c77-959a3eba0197
PlutoTest.@test std(gaussian_noise_remove(arr_rand, kernel_size=(8,8), σ=2)) < 0.05

# ╔═╡ 911c2e1e-f04a-4efa-9aff-7f4a424011ae
PlutoTest.@test sum(abs2, img_g .- img) > sum(abs2, img_g_gauss .- img)

# ╔═╡ 0f6ff33f-e50a-4317-a240-3b2f52de7a97
md"## 2.2 Noise Removal with Median Filter
The median filter slides with a quadratic shaped box over the array and
always takes the median value of this array.

So conceptually, you can use two for loops over the first two dimensions of the array. Then extract a region around this point and calculate the median.
Try to preserve the image output size. Assign the median to this pixel.
"

# ╔═╡ 7bc0f845-d5d4-4192-be1d-b97750e95761
function median_noise_remove!(arr; kernel_size=(5,5))
	# todo
	return similar(arr)
end

# ╔═╡ ac6e7ea7-c277-4180-a809-a37531ab4f11
md"kernel\_size\_2 = $(@bind ks_2 Slider(3:2:9, show_value=true))"

# ╔═╡ 8208b936-ad5d-45a0-96d1-04fd131f1e29
begin
	img_p_median = median_noise_remove!(copy(img_p), kernel_size=(ks_2,ks_2));
	img_g_median = median_noise_remove!(copy(img_g), kernel_size=(ks_2,ks_2));
	img_h_median = median_noise_remove!(copy(img_h), kernel_size=(ks_2,ks_2));
end;

# ╔═╡ 06af0b94-f4e1-4114-b565-730b87567995
my_show([img_p_median img_g_median img_h_median])

# ╔═╡ 707a3be6-ac83-421d-a8e1-31d75444aeef
md"### 2.2 Test"

# ╔═╡ 07f5e6da-7868-4e22-90aa-7a81c2b41ab3
PlutoTest.@test sum(abs2, img_p .- img) > sum(abs2, img_p_median .- img)

# ╔═╡ 3dbe2b38-8609-41e5-b3b9-d62da394b89d
PlutoTest.@test sum(abs2, img_g .- img) > sum(abs2, img_g_median .- img)

# ╔═╡ 867d5676-8dfc-4bb5-b555-73aa767f7e9e
md"## 3 Final Images"

# ╔═╡ 3605fda7-078e-4bff-9e33-91c4025e775a
my_show([img_g img_p img_h])

# ╔═╡ c09b3470-95bd-4f49-b64b-26459c95ab19
my_show([median_noise_remove!(copy(img_g)) median_noise_remove!(copy(img_p)) median_noise_remove!(copy(img_h))])  

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
ImageFiltering = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
ImageShow = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
IndexFunArrays = "613c443e-d742-454e-bfc6-1d7f8dd76566"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PoissonRandom = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
TestImages = "5e47fb64-e119-507b-a336-dd2b206d9990"

[compat]
Colors = "~0.12.8"
ImageFiltering = "~0.7.2"
ImageShow = "~0.3.6"
IndexFunArrays = "~0.2.5"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.44"
PoissonRandom = "~0.4.1"
TestImages = "~1.7.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "5260083a9dfaa0285cb23a3321980040be762f9f"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

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

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["ArrayInterfaceCore", "Compat", "IfElse", "LinearAlgebra", "Static"]
git-tree-sha1 = "d6173480145eb632d6571c148d94b9d3d773820e"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "6.0.23"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e9f7992287edfc27b3cbe0046c544bace004ca5b"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.22"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "1dd4d9f5beebac0c03446918741b1a03dc5e5788"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.6"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

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
git-tree-sha1 = "3ca828fe1b75fa84b021a7860bd039eaea84d2f2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.3.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

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
git-tree-sha1 = "5158c2b41018c5f7eb1470d558127ac274eca0c9"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.1"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

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

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78e2c69783c9753a91cdae88a8d432be85a2ab5e"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

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

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

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

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "8b251ec0582187eff1ee5c0220501ef30a59d2f7"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.2"

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

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

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

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

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

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

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

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

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

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

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

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "6e33d318cf8843dade925e35162992145b4eb12f"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.44"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "9ac1bb7c15c39620685a3a7babc0651f5c64c35b"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

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

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "de4f0a4f049a4c87e4948c04acff37baf1be01a6"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.7.7"

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
git-tree-sha1 = "3cbfd92ae1688129914450ff962acfc9ced42520"
uuid = "5e47fb64-e119-507b-a336-dd2b206d9990"
version = "1.7.0"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "70e6d2da9210371c927176cb7a56d41ef1260db7"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.1"

[[deps.TiledIteration]]
deps = ["ArrayInterface", "OffsetArrays"]
git-tree-sha1 = "5e02b75701f1905e55e44fc788bd13caedb5a6e3"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.4.1"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

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

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─31cc54e9-d47b-4df7-a785-dad0120e8cdb
# ╠═8a11c2b7-b6e5-4c26-9e06-ad3ab958db9c
# ╠═50044e20-1dda-444f-9011-d9a55e182e11
# ╟─b731d80c-3012-49f5-8e63-9b7cea0ed667
# ╟─d8c83f2c-76ea-42b9-bd39-7fa5cafcd4d3
# ╠═334ebdc5-dc11-4e83-ad12-a1961420a97b
# ╟─2122c12f-0832-4f19-9a77-3047d5311cba
# ╟─a7d7a415-7464-43fb-a2bf-8e3f2d999c32
# ╠═de104b0a-28a4-42de-ae8d-fa857a8f32e0
# ╟─f89bc1e4-d406-4550-988b-71496b64035a
# ╠═c0b5ff08-7342-4c25-9791-be3a4918c400
# ╟─d56bc500-f5cf-40f3-8c0b-b4ae26a58367
# ╠═1116455e-4698-489a-add2-4e1d8ecbb546
# ╟─cc7bbb55-b8dd-454b-b170-45137ab0c0b5
# ╠═c997c248-fba8-4ed8-ab75-179d248e6b83
# ╠═07bc00b7-81ad-4645-8fbe-8545ccdfe7e4
# ╠═d317ed1a-1867-4293-8db4-51352212929c
# ╟─e69a8b77-833a-4498-8bc9-e4c970f6529d
# ╠═3ccbf4fb-0707-4ce8-83d8-b1af555f3fe1
# ╠═1e4b32cb-00ea-439f-900e-abadc850be95
# ╠═6cf8ff2a-4255-4c53-9376-5a0d2d372569
# ╠═9a3feb13-0424-448c-934b-45c476cb0096
# ╟─e3fcd381-9652-4e1d-8585-a91f87b73ce5
# ╠═f3efcfcf-4516-43d2-9375-d8863b4b6f5b
# ╠═5b68030f-557e-449e-8baf-cf7ae65e4425
# ╠═78a39a77-bc9d-427a-a113-4e64b0c9d311
# ╠═41118b3c-3e49-4fc7-8c80-7906b7cf91fe
# ╟─f8be3a4b-d8f0-40cb-869e-81c17e58327a
# ╠═65a8981d-890b-4590-a913-6890ecf3817d
# ╠═272c8afb-667f-49c7-af48-2feb1b9d06f7
# ╟─59482e49-d784-4303-9487-bf37f6a0462e
# ╠═6e050fc2-4db2-4e6e-abd4-2812b47c070f
# ╠═22985061-9740-45e9-af56-e0787dadba7b
# ╠═b6c25124-3165-4c85-9b1c-da68bd68b59f
# ╠═76e83c0b-e9a8-46bd-8a7e-cac39ccaefff
# ╟─42f92a8c-54bb-4033-9ef2-f5ffb78f4f3a
# ╠═a9ad9d52-53c9-44ce-8edc-af8d7bfdc81f
# ╟─373a3e39-fd8c-4cea-88bb-2b9c3734f14f
# ╟─2c274829-ef33-4f34-9ba0-63d3b8cd7343
# ╠═52ebebfc-940b-47e2-a059-40ec996f2f98
# ╠═2d6f25a5-b902-4e46-ac3f-cf452f525912
# ╠═62c6e2f0-0007-40ec-81c5-61cda80f59e5
# ╠═75aa8212-8ba3-4fba-b8b2-af2dd9c82adf
# ╠═ded0b153-3ba6-4ad6-906b-87930e2b82d0
# ╠═6acf7ba6-445d-4ac0-a83f-6916318c783c
# ╟─4f683117-9a09-4837-b0af-5818ac2e62fd
# ╟─01203c2c-4e67-492c-9d9c-dc84a84afe15
# ╠═daae9a7d-99f4-4d34-a7fd-886fa754d59b
# ╠═92d43328-29a2-45b2-ae71-e169dec02316
# ╟─512486eb-1919-446a-9f24-e9713546d237
# ╠═9b83d92d-f838-4329-b0b1-62bf7a18f123
# ╠═bbbed1b7-a406-44a7-8c77-959a3eba0197
# ╠═911c2e1e-f04a-4efa-9aff-7f4a424011ae
# ╟─0f6ff33f-e50a-4317-a240-3b2f52de7a97
# ╠═7bc0f845-d5d4-4192-be1d-b97750e95761
# ╟─ac6e7ea7-c277-4180-a809-a37531ab4f11
# ╠═8208b936-ad5d-45a0-96d1-04fd131f1e29
# ╠═06af0b94-f4e1-4114-b565-730b87567995
# ╟─707a3be6-ac83-421d-a8e1-31d75444aeef
# ╠═07f5e6da-7868-4e22-90aa-7a81c2b41ab3
# ╠═3dbe2b38-8609-41e5-b3b9-d62da394b89d
# ╟─867d5676-8dfc-4bb5-b555-73aa767f7e9e
# ╠═3605fda7-078e-4bff-9e33-91c4025e775a
# ╠═c09b3470-95bd-4f49-b64b-26459c95ab19
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
