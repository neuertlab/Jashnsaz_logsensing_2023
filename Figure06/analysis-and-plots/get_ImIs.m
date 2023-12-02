function [time,Im,Is] = get_ImIs(data,tmin,tmax,type)
    index = find(data.tt>=tmin & data.tt<=tmax); N = numel(index); 
    if type=="Vol"
        time = data.tt(index); Im = data.Volm(index); Is = data.Vols(index);
    elseif type=="Hog1"
        try
            time = data.tt(index); Im = data.mHog1nuc(index); Is = data.sHog1nuc(index);
        catch
            time = data.tt(index); Im = data.Hog1m(index); Is = data.Hog1s(index);
        end
    end
    time=reshape(time,N,1); Im=reshape(Im,N,1); Is=reshape(Is,N,1); 
end