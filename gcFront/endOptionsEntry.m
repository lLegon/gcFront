function endOptionsEntry(input1,~)
% function used to resume gcFront after option entry graphic is closed

if isequal(class(input1),'matlab.ui.control.Button')
    % function was called by pressing the start algorithm button- get
    % parent figure
    input1=input1.Parent;
end



uiresume(input1);




end