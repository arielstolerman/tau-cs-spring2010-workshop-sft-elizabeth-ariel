% runs the SFT algorithm (the sft_dp.m version) for numOfRec times and returns the intersection of the
% results. Using this script instead of running only one reccurence of the SFT clears noise and false-significant elements.
% Returns:
% - L, coeffs - as described in the original SFT script, only this is an intersection of results from numOfRec reccurences.
% - sizes - a vector of the size of L (and coeffs) in each reccurence.

function[L,coeffs,sizes]=sft_dp_rec(isLogged,G,tau,func,m_A,m_B,numOfIterations,numOfRec)

% hold the size of each result
sizes = zeros(1,numOfRec);

[L,coeffs]=sft_dp(isLogged,G,tau,func,m_A,m_B,numOfIterations);
sizes(1) = size(L,1);
% run numOfRec - 1 more times
for i=2:numOfRec
    [tmpL,tmpCoeffs]=sft_dp(isLogged,G,tau,func,m_A,m_B,numOfIterations);
    indToRemove = []; % vector to hold indices of elements to be removed from L
    ind = 1;
    for j=1:size(L,1)
        removeElem = 1; % remove current element in L unless found a match
        for k=1:size(tmpL,1)
            if (min(L(j,:)==tmpL(k,:)))
                % found match, do NOT remove the element from L
                removeElem = 0;
                break;
            end
        end
        % mark to remove the element from L
        if (removeElem)
            indToRemove(ind) = j;
            ind = ind + 1;
        end
    end
    % remove all marked indices in L and coeffs
    indToRemove = sort(indToRemove,'descend');
    for j=1:length(indToRemove)
        L(indToRemove(j),:) = [];
        coeffs(indToRemove(j),:) = [];
    end
    % save current iteration size of L
    sizes(i) = size(L,1);
end