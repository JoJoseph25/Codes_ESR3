function checkBoxCallback(source, data, index)
%     fprintf('Check box %s is %d\n',get(source,'String'),get(source,'Value'));
    check=get(source,'Value');
    if check==1
        set(handles.checkbox2,'Visible','off')
    else
        set(handles.checkbox2,'Visible','on')
    end
end
