function aaCrossCorr(x, y, n, NumLags)
    arguments
        x
        y
        n
        NumLags = 40
    end

    figure;
    stem(-NumLags:NumLags, crosscorr(x, y, NumLags=NumLags));
    title('Crosscorrelation function');
    xlabel('Lag');
    
    hold on
    plot(-NumLags:NumLags, 2/sqrt(n) * ones(1, 2*NumLags+1), '--')
    plot(-NumLags:NumLags, -2/sqrt(n) * ones(1, 2*NumLags+1), '--')
    hold off
end