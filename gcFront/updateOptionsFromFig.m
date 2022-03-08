function updateOptionsFromFig(textEntry,~,field)
% during option entry by user, read the value that is in the options text 
% box, and store this value in the options struct

if isempty(field)
    return
end

% get option figure and value to store
optionFig=textEntry.Parent;
entryValue=join( string(textEntry.Value), ''); %#ok<*NASGU>

% update the user data in the options figure with whatever has just been
% entered into the text box
eval(['optionFig.UserData.',field,'=entryValue;']);

end