function ruleMat=readRulesField(testRule)
% reads a generule from a metabolic model and converts it into a matrix of
% gene indices, where columns represent a set of genes that can carry out
% the reaction. 

origTestRule=testRule;

% remove spaces from the rule
testRule=strrep(testRule,' ','');

% return if not gene associated
if isempty(testRule)
    ruleMat=[];
    return
end

% stop if rule contains unpaired bracket, as rule won't work otherwise
if count(testRule,'(')~=count(testRule,')')
    ruleMat=[];
    warning(['Rules field has an unpaired brackets: ',origTestRule])
    return
end

% replace the format 'x(geneNumber)' with just the gene number
testRule=strrep(testRule,'x(','');
numberInds=regexp(testRule,'\d');
toDel=false(size(testRule));
for a=1:length(numberInds)
    if a==length(numberInds) || numberInds(a)+1~=numberInds(a+1)
        % the next character is a close bracket that needs to be removed
        toDel(numberInds(a)+1)=true;
    end
end
testRule(toDel)='';


bracketMap=containers.Map;
bCount=1;

% break down any brackets into a matrix
while 1
    % find the indices of the text to analyse (the contents of a bracket pair, or if no
    % brackets, then entire text)
    startBrack=0;
    endBrack=length(testRule)+1;
    for a=1:length(testRule)
        if testRule(a)=='('
            startBrack=a;
        elseif testRule(a)==')'
            endBrack=a;
            break
        end
    end
    
    
    % get the text to be analysed
    containedText=testRule(startBrack+1:endBrack-1);
    
    if contains(containedText,'|')
        % combining genes with an 'OR' statement
        
        % check that rule is not ambiguous
        if contains(containedText,'&')
            % text contains both and + or in the same statement, with no
            % brackets to clarify which operation to carry out first
            warning(['Rules are ambiguous, they contain | and & without a bracket to specify what they refer to:', origTestRule])
            ruleMat=[];
            return
        end
        
        % split text along "or"
        % must replace the | character as splitString doesn't like
        % splitting this character for some reason
        containedText=strrep(containedText,'|','*');
        splitText=splitString(containedText,'*')';
        
        % take all things that aren't a bracket, and put them into a matrix
        % of gene indices
        notAt=~contains(splitText,'@');
        atInds=find(~notAt);
        
        bracketMatrix=str2double(string(splitText(notAt)));
        
        % add contents of brackets as extra columns of matrix
        if ~isempty(atInds)
            
            for a=1:length(atInds)
                
                newInds=bracketMap(splitText{atInds(a)});
                
                % add infinity into spaces in case there is a size mismatch between
                % brackets and indices matrix
                if size(bracketMatrix,1)<size(newInds,1)
                    bracketMatrix=[bracketMatrix;ones(size(newInds,1)-size(bracketMatrix,1), size(bracketMatrix,2))*inf];
                elseif size(bracketMatrix,1)>size(newInds,1)
                    newInds=[newInds; ones(size(bracketMatrix,1)-size(newInds,1), size(newInds,2))*inf];
                end
                
                bracketMatrix=[bracketMatrix,newInds];
                
            end
        end
        

    elseif contains(containedText,'&')
        % Combining genes with an 'AND' statement
        
        splitText=splitString(containedText,'&');
        
        % take all things that aren't a bracket, and put them into a matrix
        % of gene indices
        notAt=~contains(splitText,'@');
        atInds=find(~notAt);
        
        bracketMatrix=str2double(string(splitText(notAt)));
        
        % add contents of brackets as extra rows in matrix. If bracket
        % contained an 'OR' statement, this will create all possible
        % combinations of bracket contents and things outside of it
        if ~isempty(atInds)
            
            
            for a=1:length(atInds)
                
                bracketMatrix=combvec(bracketMatrix,bracketMap(splitText{atInds(a)}));
                
            end
        end
        
    else
        % no 'AND' or 'OR' statements in brackets
        
        if contains(containedText,'@')
            bracketMatrix=bracketMap(containedText);
        else
            bracketMatrix=str2double(string(containedText));
        end
        
        
    end
    

    % assign a new name that refers to the analysed text
    newName=['@',num2str(bCount)];
    
    % store gene indices that this name refers to in map
    % container
    bracketMap(newName)=bracketMatrix;
    bCount=bCount+1;
    
    % replace analysed text with its new name
    if startBrack~=0
        testRule=strrep(testRule,testRule(startBrack:endBrack),newName);
    else
        testRule=newName;
    end
    
    % if all logical statements have been analysed, terminate
    if ~contains(testRule,'|') && ~contains(testRule,'&')
        break
    end
end

% retrieve the matrix that corresponds to the full text
ruleMat=bracketMap(testRule);

% sort matrix and remove non-unique columns
ruleMat=sort(ruleMat);
ruleMat=unique(ruleMat','rows')';


end