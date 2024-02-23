# Repeat functions to repeat patterns 

function repeat_pattern(a::Vector{Int64}) 
    
    # Find the length of the last non-negative segment
    last_value_length = 0
    for i in length(a):-1:1
        if a[i] != -1
            last_value_length = length(a) - i
            break
        end
    end

	# Insert the right switch in each swicth setup i.e, fill out the holes up to that point
	for i in 1:last_value_length-1
		if a[i] == -1
			a[i] = a[i-1]
		end
	end

	# find the length of the second to last switch 
	# The length of the last switch will be set the same length 
	length_last_switch = count(x -> x == a[last_value_length-2], a)
	a[last_value_length-1:last_value_length-2+length_last_switch] .= a[last_value_length-1]
    
	# now it should repeat the pattern as many times as possible 
	idx = last_value_length - 2 + length_last_switch # new last non-zero value 

	nrows = size(a)[1]
	pattern = view(a, 1:idx)
	repetitions = div(nrows, idx)
	remainder = mod(nrows, idx)
	
	repeated_pattern = repeat(pattern, repetitions, 1)
	remainder_pattern = view(pattern, 1:remainder, :)
	
	return vcat(repeated_pattern, remainder_pattern)
end

# Setting up the switches such that they follow a repeated pattern. 
# that means if the section times are longer than the number of switches, the switches should be repeated 
# such that one does not need to specify all the switches in for example in an SMB but only one cycle. 
# it is done so in a row-wise manner such that if a switch is not determined at a specific section time, it will be overwritten by the cycle
# The pattern is repeated of vector A and onto vector B itself. 
function repeat_pattern!(B::Vector, A::Array)
    idx = findfirst(x -> x == -1, A[:, 1])
    if idx !== nothing
        nrows = size(B)[1]
        
        # Extract the pattern to repeat
        pattern = B[1:idx-1]
        
        # Compute how many times the pattern can fit in the remaining part of the matrix
        pattern_repetitions = div(nrows - (idx-1), size(pattern, 1))
        
        # Repeat the pattern as many times as possible
        repeated_pattern = repeat(pattern, pattern_repetitions, 1)
        
        # Calculate the remaining rows
        remaining_rows = nrows - (idx-1) - size(repeated_pattern, 1)
        
        # If there are remaining rows, append the pattern again as needed
        if remaining_rows > 0
            extra_rows = repeat(pattern[1:remaining_rows, :], 1)
            repeated_pattern = vcat(repeated_pattern, extra_rows)
        end
        
        # Replace the content of B with the repeated pattern
        B .= vcat(B[1:idx-1, :], repeated_pattern)
    end
    return B
end
function repeat_pattern!(B::Matrix, A::Array)
    idx = findfirst(x -> x == -1, A[:, 1])
    if idx !== nothing
        nrows, ncols = size(B)
        
        # Extract the pattern to repeat
        pattern = B[1:idx-1, :]
        
        # Compute how many times the pattern can fit in the remaining part of the matrix
        pattern_repetitions = div(nrows - (idx-1), size(pattern, 1))
        
        # Repeat the pattern as many times as possible
        repeated_pattern = repeat(pattern, pattern_repetitions, 1)
        
        # Calculate the remaining rows
        remaining_rows = nrows - (idx-1) - size(repeated_pattern, 1)
        
        # If there are remaining rows, append the pattern again as needed
        if remaining_rows > 0
            extra_rows = repeat(pattern[1:remaining_rows, :], div(ncols, size(pattern, 2)), 1)
            repeated_pattern = vcat(repeated_pattern, extra_rows)
        end
        
        # Replace the content of B with the repeated pattern
        B .= vcat(B[1:idx-1, :], repeated_pattern)
    end
    return B
end

function repeat_pattern!(B::Array, A::Array)
    idx = findfirst(x -> x == -1, A[:, 1])
    if idx !== nothing
        nrows, ncols, _ = size(B)
        
        # Extract the pattern to repeat
        pattern = B[1:idx-1, :, :]
        
        # Compute how many times the pattern can fit in the remaining part of the matrix
        pattern_repetitions = div(nrows - (idx-1), size(pattern, 1))
        
        # Repeat the pattern as many times as possible
        repeated_pattern = repeat(pattern, pattern_repetitions, 1)
        
        # Calculate the remaining rows
        remaining_rows = nrows - (idx-1) - size(repeated_pattern, 1)
        
        # If there are remaining rows, append the pattern again as needed
        if remaining_rows > 0
            extra_rows = repeat(pattern[1:remaining_rows, :, :], div(ncols, size(pattern, 2)), 1)
            repeated_pattern = vcat(repeated_pattern, extra_rows)
        end
        
        # Replace the content of B with the repeated pattern
        B .= vcat(B[1:idx-1, :, :], repeated_pattern)
    end
    return B
end

function repeat_pattern!(B::Matrix, idx::Int64)
    nrows, ncols = size(B)
    
    # Extract the pattern to repeat
    pattern = B[1:idx, :]
    
    # Compute how many times the pattern can fit in the remaining part of the matrix
    pattern_repetitions = div(nrows - idx, size(pattern, 1))
    
    # Repeat the pattern as many times as possible
    repeated_pattern = repeat(pattern, pattern_repetitions, 1)
    
    # Calculate the remaining rows
    remaining_rows = nrows - idx - size(repeated_pattern, 1)
    
    # If there are remaining rows, append the pattern again as needed
    if remaining_rows > 0
        extra_rows = repeat(pattern[1:remaining_rows, :], div(ncols, size(pattern, 2)), 1)
        repeated_pattern = vcat(repeated_pattern, extra_rows)
    end
    
    # Replace the content of B with the repeated pattern
    B .= vcat(B[1:idx, :], repeated_pattern)
end




