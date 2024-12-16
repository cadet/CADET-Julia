# Repeat functions to repeat patterns 


function repeat_pattern(a::Vector{Int64}) 
    
    # Find the length of the last non-negative segment
    last_value_length = 0
    for i in length(a):-1:1
        if a[i] != -1
            last_value_length = i
            break
        end
    end

	# Insert the right switch in each swicth setup i.e, fill out the holes up to that point
	for i in 1:last_value_length
		if a[i] == -1
			a[i] = a[i-1]
		end
	end

	# find the length of the second to last switch 
	# The length of the last switch will be set the same length 
	length_last_switch = count(x -> x == a[last_value_length-1], a)
	a[last_value_length:last_value_length-1+length_last_switch] .= a[last_value_length]
    
	# now it should repeat the pattern as many times as possible 
	idx = last_value_length - 1 + length_last_switch # new last non-zero value 

	nrows = size(a)[1]
	pattern = view(a, 1:idx)
	repetitions = div(nrows, idx)
	remainder = mod(nrows, idx)
	
	repeated_pattern = repeat(pattern, repetitions, 1)
	remainder_pattern = view(pattern, 1:remainder, :)
	
	return vec(vcat(repeated_pattern, remainder_pattern))
end

# Repeat pattern from a certain idx specification cIn = switches.ConnectionInstance.cIn_c
function repeat_pattern(cIn, switchSetup, nSwitches) 	
    
    length_last_cycle = 1
    idx = 0
    for i in eachindex(switchSetup)
        if switchSetup[i] == nSwitches
            for j = i+1:length(switchSetup)
                if switchSetup[j] == switchSetup[i]
                    length_last_cycle += 1
                else 
                    idx = i
                    break 
                end
            end
			break
        end 
    end
    
	# now it should repeat the pattern as many times as possible 
	nrows = size(cIn)[1]
	pattern = cIn[1:idx][:][:] # view(cIn, 1:idx, :, :)
	if idx == 0
		return cIn
	else
		repetitions = div(nrows, idx)
		remainder = mod(nrows, idx)

		repeated_pattern = repeat(pattern, repetitions)
		remainder_pattern = pattern[1:remainder][:][:] # view(pattern, 1:remainder, :, :)

		return vcat(repeated_pattern, remainder_pattern)
	end
end



