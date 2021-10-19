

%Plot the insects ID regarding the frames where they appear
function frames_with_insect(attr_frame, attr_id)
    %PLOT OBJ_ID OVER FRAME
    figure()
    plot(attr_frame, attr_id, '+');
    title('Relation insect ID and frame where it appears');
    xlabel('Frame Number');
    ylabel('Insect ID');
end