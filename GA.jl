### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 7f439040-12f0-11ef-35d9-630d8b2d34f3
begin
	using Random
	using HTTP
	using PlutoUI
end

# ╔═╡ f59b7c86-c2b2-444a-8917-edf8318875d3
begin
	PlutoUI.HTML("<img src=https://editor.analyticsvidhya.com/uploads/779521623127686493.png>")
end

# ╔═╡ d44c1910-6f08-44ed-ba51-5241d5ec2603
function firstGeneration(
	 numberOfChromosomes::Int
	,chromosomeLength::Int=32
	,population_size::Int=20
	,encoding::String="binary")
	if encoding == "binary"
		return [rand(0:1,numberOfChromosomes*chromosomeLength)
			for _ in 1:population_size]
	end
end

# ╔═╡ ab9ff0a5-9faa-4a86-b439-cc86fdf8aec6
# selection
function select_parents(fitnessArray::Vector,objective::String="max")
	parents::Vector{Int8} = rand(1:length(fitnessArray),2)
	fitnessArrayCopy = copy(fitnessArray)
	for i in range(1,2)
		probability = rand(Float64,1)[1]
		max_fitness = maximum(fitnessArrayCopy)
		fitnessArrayCopy = (fitnessArrayCopy ./ max_fitness)
		min_fitness = minimum(fitnessArrayCopy)
		fitnessArrayCopy = (fitnessArrayCopy .- min_fitness)
		totalFitnesses = sum(fitnessArrayCopy)
		prevPropabilityToSelect = 0.0
		propabilityToSelect = 0.0
		for (nationalId,fitness) in enumerate(fitnessArrayCopy)
			propabilityToSelect = fitness / totalFitnesses
			propabilityToSelect += prevPropabilityToSelect

			if (probability >= prevPropabilityToSelect && 
				probability <= propabilityToSelect)
				parents[i] = nationalId
				fitnessArrayCopy = deleteat!(fitnessArrayCopy,nationalId)
				break
			end
			prevPropabilityToSelect = propabilityToSelect
		end
	end
	return parents
end

# ╔═╡ 2b0ed86d-75e5-455c-b9ac-515b60161ebd
PlutoUI.HTML("<img src=https://media.geeksforgeeks.org/wp-content/uploads/20190620121313/twopointCrossover-2.png>")

# ╔═╡ 96ccbad2-ecca-4ee3-b3e5-2bcec9652b8e
function mate(father::Vector{Int}, mother::Vector{Int})
    # Ensure parents have the same length
    @assert length(father) == length(mother)

    # Select a random crossover point
    crossover_point = rand(1:length(father))

    # Perform crossover
    child1 = vcat(father[1:crossover_point], mother[crossover_point+1:end])
    child2 = vcat(mother[1:crossover_point], father[crossover_point+1:end])

    return child1, child2
end


# ╔═╡ cd2deab6-2c94-48a0-be08-12249cb50dd1
PlutoUI.HTML("<img src=https://thatgamesguy.co.uk/wp-content/uploads/2018/05/ga_crossover_mutation.png>")

# ╔═╡ 85cf9276-1c8f-4a71-9485-9e9e57057366
function mutate(child::Vector{Int}, mutation_rate::Float64)
    mutated_child = copy(child)
    for i in eachindex(mutated_child)
        if rand() < mutation_rate
            mutated_child[i] = 1 - mutated_child[i]  # Flip the bit
        end
    end
    return mutated_child
end

# ╔═╡ 465370ed-1dc8-466c-b770-e436767a2170
# Function to convert a 32-bit binary array to a decimal number using IEEE 754 standard from chatGPT
function ieee754_to_float32(binary_array::Vector{Int})
	if length(binary_array) != 32
		throw(ArgumentError("The binary array must have exactly 32 elements."))
	end

	# Extract the sign, exponent, and mantissa
	sign_bit = binary_array[1]
	exponent_bits = binary_array[2:9]
	mantissa_bits = binary_array[10:end]

	# Calculate the sign
	sign = (-1.0)^sign_bit

	# Calculate the exponent
	exponent = sum(exponent_bits[i] * 2^(8-i) for i in 1:8) - 127

	# Calculate the mantissa
	mantissa = 1.0 + sum(mantissa_bits[i] * 2.0^(-i) for i in 1:23)

	# Calculate the final float value
	float_value = sign * mantissa * 2.0^exponent
	println(float_value)
	if !isfinite(float_value)
		if isnan(float_value)
			float_value = 0
		else
			float_value = floatmax(Float32)-1
		end
	end

	return float_value
end

# ╔═╡ 67599e83-5db3-4716-9a49-3c70c2b8af40
# Function to convert a binary array to a number
function binary_to_decimal(binary_array::Vector{Int})
    decimal_number = 0
    for (i, bit) in enumerate(reverse(binary_array))
        decimal_number += bit * (2^(i-1))
    end
    return decimal_number
end

# ╔═╡ db88d4dd-2363-43ab-9966-d80b051e2c44
# Function to segment binary vector, convert segments, and evaluate function f
# begin
function calculateFitness(generation::Vector{Vector{Int}},f,
	chromosomeLength::Int=32,decoding=binary_to_decimal)
	
	fitnessArray::Vector = []
	numberOfChromosomes = div(length(generation[1]), chromosomeLength)
	for human in generation
		chromosomes = [human[(i-1)*chromosomeLength + 1:i*chromosomeLength] for i in 1:numberOfChromosomes]
		
		# Convert each chromosome to a float using decode
		chromosomes_values = [decoding(chromosome) for chromosome in chromosomes]
		# Evaluate the function f with the float values
		fitness = f(chromosomes_values...)
		push!(fitnessArray, fitness)
	end

	return fitnessArray
end

# ╔═╡ eac35b19-5412-4f22-8f2d-d77b56e0b108
function genetic_algorithm(
	 numberOfChromosomes::Int
	,f
	,population_size::Int=20
	,chromosomeLength::Int=32
	,mutation_rate::Float64 = 0.001
	,encoding::String="binary"
	,decoding=binary_to_decimal
	,stopAt::Int64 = 10000
	,objective::String = "max"
)
	generation = firstGeneration(2,chromosomeLength,20)
	fitnesses = calculateFitness(generation,f,chromosomeLength)
	for i in range(1,stopAt)
		fatherId, motherId = select_parents(fitnesses,objective)
		father = generation[fatherId]
		mother = generation[motherId]
		child1, child2 = mate(father,mother)
		child1 = mutate(child1,mutation_rate)
		child2 = mutate(child2,mutation_rate)
		victimId = argmin(fitnesses)
		generation[victimId] = child1
		victimId = argmin(fitnesses)
		generation[victimId] = child2
		fitnesses = calculateFitness(generation,f,chromosomeLength)
	end
	fittestId = argmax(fitnesses)
	fittest = generation[fittestId]
	chromosomes = [fittest[(i-1)*chromosomeLength + 1:i*chromosomeLength] for i in 1:numberOfChromosomes]
	fittest_chromosome = [decoding(chromosome) for chromosome in chromosomes]
	return fittest_chromosome
end


# ╔═╡ 46c7eb2f-cc24-430f-9731-cf3f66bd5bb5
begin
	# 	function genetic_algorithm(
	# 	 numberOfChromosomes::Int
	# 	,f
	# 	,population_size::Int=20
	# 	,chromosomeLength::Int=32
	# 	,mutation_rate::Float64 = 0.001
	# 	,encoding::String="binary"
	# 	,decoding=binary_to_decimal
	# 	,stopAt::Int64 = 10000
	# 	,objective::String = "max"
	# )
	# Note number of chormosomes is the number of variables
	f(x,y) = sin(x^2 + y^2) * cos(x * y) - log(x^2 + y^2 + 1)
	genetic_algorithm(
    2,
    f,
    20,
    31,
    0.005,
    "binary",
    binary_to_decimal,
    20000,
    "min"
    )
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
HTTP = "~1.10.8"
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "006f0907a31010c30126de19f6f7be2568e3dd35"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "6cbbd4d241d7e6579ab354737f4dd95ca43946e1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3da7367955dcc5c54c1ba4d402ccdc09a1a3e046"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "5d54d076465da49d6746c647022f3b3674e64156"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.8"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═7f439040-12f0-11ef-35d9-630d8b2d34f3
# ╠═f59b7c86-c2b2-444a-8917-edf8318875d3
# ╠═d44c1910-6f08-44ed-ba51-5241d5ec2603
# ╠═db88d4dd-2363-43ab-9966-d80b051e2c44
# ╠═ab9ff0a5-9faa-4a86-b439-cc86fdf8aec6
# ╠═2b0ed86d-75e5-455c-b9ac-515b60161ebd
# ╠═96ccbad2-ecca-4ee3-b3e5-2bcec9652b8e
# ╠═cd2deab6-2c94-48a0-be08-12249cb50dd1
# ╠═85cf9276-1c8f-4a71-9485-9e9e57057366
# ╠═eac35b19-5412-4f22-8f2d-d77b56e0b108
# ╠═46c7eb2f-cc24-430f-9731-cf3f66bd5bb5
# ╠═465370ed-1dc8-466c-b770-e436767a2170
# ╠═67599e83-5db3-4716-9a49-3c70c2b8af40
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
