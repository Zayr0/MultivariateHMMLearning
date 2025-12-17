function annotatedSleepStages = TransformAnnotatedData(Annotations, blockTime)
    AnnotatedStages = Annotations.Annotations.Annotations;
    AnnotatedStages(AnnotatedStages == "Sleep stage W") = 1;
    AnnotatedStages(AnnotatedStages == "Sleep stage 1") = 2;
    AnnotatedStages(AnnotatedStages == "Sleep stage 2") = 3;
    AnnotatedStages(AnnotatedStages == "Sleep stage 3") = 4;
    AnnotatedStages(AnnotatedStages == "Sleep stage 4") = 4;
    AnnotatedStages(AnnotatedStages == "Sleep stage R") = 5;

    AnnotatedStages = double(AnnotatedStages);
    
    Durations = seconds(Annotations.Annotations.Duration);
    
    annotatedSleepStages = [];
    for i = 1:length(Durations)-1
        annotatedSleepStages = [annotatedSleepStages; AnnotatedStages(i) * ones(round(Durations(i)/blockTime), 1)];
    end
end