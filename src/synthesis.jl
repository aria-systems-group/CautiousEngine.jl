
function general_label_fcn(state_extent, default_label, unsafe_label, labels_dict; unsafe=false)
    if unsafe
        return unsafe_label 
    end
    state_label = default_label
    for label in keys(labels_dict) 
        for extent in labels_dict[label]
            flags = [sum(state_extent[dim])/2∈extent[dim] for dim in keys(extent)]
            if sum(flags) == length(keys(extent))
                state_label = label
                @debug state_label
                break
            end
        end
    end
    return state_label
end
