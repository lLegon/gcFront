function printYieldoverMu_add(muRange,yieldRange,prodRange,metName)
    grid on
    
    title(metName);
    
    [hAx1,hLine11,hLine21]     = plotyy(muRange,yieldRange,muRange,prodRange);

    hLine21.LineStyle   = '--';   
    hLine21.Color       = 'black'; 
    hLine11.Color       = 'black';    

    set(hAx1,'YColor','black');  
    set(get(hAx1(1),'Ylabel'),'String','Production rate [mmol/g/h]');
    set(get(hAx1(2),'Ylabel'),'String','Productivity [mol/h]');

    xlabel('Growth Rate [1/h]');    ylabel('Production rate [mmol/g/h]');
end