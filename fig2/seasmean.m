function seasonal_mean = seasmean(x,season)
% x needs to be an array that starts with January and must be a multiple of 12

if (size(season) == [1 4])

    if (season == 'DJFM')
        for j=2:length(x)/12
            seasonal_mean(j) = ((31*x(12*(j-1)+12))+(31*x(12*(j-1)+1))+(28*x(12*(j-1)+2))+(31*x(12*(j-1)+3)))/121 ;
        end
        seasonal_mean(1) = NaN;
        return
    end

    if (season == 'MAMJ')
        for j=1:length(x)/12
            seasonal_mean(j) = ((31*x(12*(j-1)+3))+(30*x(12*(j-1)+4))+(31*x(12*(j-1)+5))+(30*x(12*(j-1)+6)))/122 ;
        end
        return
    end

    if (season == 'AMJJ')
        for j=1:length(x)/12
            seasonal_mean(j) = ((30*x(12*(j-1)+4))+(31*x(12*(j-1)+5))+(30*x(12*(j-1)+6))+(31*x(12*(j-1)+7)))/122 ;
        end
        return
    end

    if (season == 'JJAS')
        for j=1:length(x)/12
            seasonal_mean(j) = ((30*x(12*(j-1)+6))+(31*x(12*(j-1)+7))+(31*x(12*(j-1)+8))+(30*x(12*(j-1)+9)))/122 ;
        end
        return
    end

    if (season == 'JFMA')
        for j=1:length(x)/12
            seasonal_mean(j) = ((31*x(12*(j-1)+1))+(28*x(12*(j-1)+2))+(31*x(12*(j-1)+3))+(30*x(12*(j-1)+4)))/120 ;
        end
        return
    end

elseif (size(season) == [1 9])

    if (season == 'SONDJFMAM')
        for j=1:(length(x)/12)-1
            seasonal_mean(j) = ((30*x(12*(j-1)+9))+(31*x(12*(j-1)+10))+(30*x(12*(j-1)+11))+(31*x(12*(j-1)+12))+(31*x(12*(j)+1))+(28*x(12*(j)+2))+(31*x(12*(j)+3))+(30*x(12*(j)+4))+(31*x(12*(j)+5)))/273 ;
        end
%        seasonal_mean(1) = seasonal_mean(2);
        return
    end

elseif (size(season) == [1 5])

    if (season == 'FMAMJ')
        for j=1:length(x)/12
            seasonal_mean(j) = ((28*x(12*(j-1)+2))+(31*x(12*(j-1)+3))+(30*x(12*(j-1)+4))+(31*x(12*(j-1)+5))+(30*x(12*(j-1)+6)))/150 ;
        end
        return
    end

    if (season == 'JFMAM')
        for j=1:length(x)/12
            seasonal_mean(j) = ((31*x(12*(j-1)+1))+(28*x(12*(j-1)+2))+(31*x(12*(j-1)+3))+(30*x(12*(j-1)+4))+(31*x(12*(j-1)+5)))/151 ;
        end
        return
    end

    if (season == 'NDJFM')
        for j=2:(length(x)/12)-1
            seasonal_mean(j) = ((30*x(12*(j-1)+11))+(31*x(12*(j-1)+12))+(31*x(12*(j)+1))+(28*x(12*(j)+2))+(31*x(12*(j)+3)))/151 ;
        end
        seasonal_mean(1) = seasonal_mean(2);
        return
    end

    if (season == 'MAMJJ')
        for j=1:(length(x)/12)
            seasonal_mean(j) = ((31*x(12*(j-1)+3))+(30*x(12*(j-1)+4))+(31*x(12*(j-1)+5))+(30*x(12*(j-1)+6))+(31*x(12*(j-1)+7)))/153 ;
        end
        seasonal_mean(1) = seasonal_mean(2);
        return
    end


elseif (size(season) == [1 6])

    if (season == 'JFMAMJ')
        for j=1:length(x)/12
            seasonal_mean(j) = ((31*x(12*(j-1)+1))+(28*x(12*(j-1)+2))+(31*x(12*(j-1)+3))+(30*x(12*(j-1)+4))+(31*x(12*(j-1)+5))+(30*x(12*(j-1)+6)))/181 ;
        end
        return
    end

    if (season == 'FMAMJJ')
        for j=1:length(x)/12
            seasonal_mean(j) = ((28*x(12*(j-1)+2))+(31*x(12*(j-1)+3))+(30*x(12*(j-1)+4))+(31*x(12*(j-1)+5))+(30*x(12*(j-1)+6))+(31*x(12*(j-1)+7)))/181 ;
        end
        return
    end

    if (season == 'MAMJJA')
        for j=1:length(x)/12
            seasonal_mean(j) = ((31*x(12*(j-1)+3))+(30*x(12*(j-1)+4))+(31*x(12*(j-1)+5))+(30*x(12*(j-1)+6))+(31*x(12*(j-1)+7))+(31*x(12*(j-1)+8)))/184 ;
        end
        return
    end

    if (season == 'NDJFMA')
      for j=2:(length(x)/12)-1
            seasonal_mean(j) = ((30*x(12*(j-1)+11))+(31*x(12*(j-1)+12))+(31*x(12*(j)+1))+(28*x(12*(j)+2))+(31*x(12*(j)+3))+(30*x(12*(j)+4)))/181 ;
        end
        seasonal_mean(1) = seasonal_mean(2);
        return
    end

    if (season == 'SONDJF')
      for j=1:(length(x)/12)-1
            seasonal_mean(j+1) = ((30*x(12*(j-1)+9))+(31*x(12*(j-1)+10))+(30*x(12*(j-1)+11))+(31*x(12*(j-1)+12))+(31*x(12*(j)+1))+(28*x(12*(j)+2)))/181 ;
        end
        seasonal_mean(1) = NaN;
        return
    end

    if (season == 'AMJJAS')
      for j=1:(length(x)/12)
            seasonal_mean(j) = ((30*x(12*(j-1)+4))+(31*x(12*(j-1)+5))+(30*x(12*(j-1)+6))+(31*x(12*(j-1)+7))+(31*x(12*(j-1)+8))+(30*x(12*(j-1)+9)))/183 ;
        end
        return
    end

    if (season == 'ONDJFM')
      for j=1:(length(x)/12)-1
            seasonal_mean(j+1) = ((31*x(12*(j-1)+10))+(30*x(12*(j-1)+11))+(31*x(12*(j-1)+12))+(31*x(12*j+1))+(28*x(12*j+2))+(31*x(12*j+3)))/182 ;
      end
      seasonal_mean(1) = NaN;
      return
    end

elseif (size(season) == [1 7])

    if (season == 'JFMAMJJ')
        for j=1:length(x)/12
            seasonal_mean(j) = ((31*x(12*(j-1)+1))+(28*x(12*(j-1)+2))+(31*x(12*(j-1)+3))+(30*x(12*(j-1)+4))+(31*x(12*(j-1)+5))+(30*x(12*(j-1)+6))+(31*x(12*(j-1)+6)))/212 ;
        end
        return
    end

elseif (size(season) == [1 3])

    if (season == 'DJF')
        for j=2:length(x)/12
            seasonal_mean(j) = ((31*x(12*(j-1)))+(31*x(12*(j-1)+1))+(28*x(12*(j-1)+2)))/90 ;
            %seasonal_mean(j) = (x(12*(j-1))+x(12*(j-1)+1)+x(12*(j-1)+2))/3 ;
        end
        seasonal_mean(1) = NaN;%seasonal_mean(2);
        return
    else

    if (season == 'JFM')
        for j=1:length(x)/12
            seasonal_mean(j) = ((31*x(12*(j-1)+1))+(28*x(12*(j-1)+2))+(31*x(12*(j-1)+3)))/90 ;
        end
        return
    else

    if (season == 'FMA')
        for j=1:length(x)/12
            seasonal_mean(j) = ((28*x(12*(j-1)+2))+(31*x(12*(j-1)+3))+(30*x(12*(j-1)+4)))/89 ;
        end
        return
    else

    if (season == 'MAM')
        for j=1:length(x)/12
            seasonal_mean(j) = ((31*x(12*(j-1)+3))+(30*x(12*(j-1)+4))+(31*x(12*(j-1)+5)))/92 ;
        end
        return
    else

    if (season == 'AMJ')
        for j=1:length(x)/12
            seasonal_mean(j) = ((30*x(12*(j-1)+4))+(31*x(12*(j-1)+5))+(30*x(12*(j-1)+6)))/91 ;
        end
        return
    else

    if (season == 'MJJ')
        for j=1:length(x)/12
            seasonal_mean(j) = ((31*x(12*(j-1)+5))+(30*x(12*(j-1)+6))+(31*x(12*(j-1)+7)))/92 ;
        end
        return
    else

    if (season == 'JJA')
        for j=1:length(x)/12
            seasonal_mean(j) = ((30*x(12*(j-1)+6))+(31*x(12*(j-1)+7))+(31*x(12*(j-1)+8)))/92 ;
        end
        return
    else

    if (season == 'JAS')
        for j=1:length(x)/12
            seasonal_mean(j) = ((31*x(12*(j-1)+7))+(31*x(12*(j-1)+8))+(30*x(12*(j-1)+9)))/92 ;
        end
        return
    else

    if (season == 'SON')
        for j=1:length(x)/12
            seasonal_mean(j) = ((30*x(12*(j-1)+9))+(31*x(12*(j-1)+10))+(30*x(12*(j-1)+11)))/91 ;
        end
        return
    else

    if (season == 'OND')
        for j=1:length(x)/12
            seasonal_mean(j) = ((31*x(12*(j-1)+9))+(30*x(12*(j-1)+10))+(31*x(12*(j-1)+11)))/92 ;
        end
        return
    else

    if (season == 'JFM')
        for j=1:length(x)/12
            seasonal_mean(j) = ((31*x(12*(j-1)+1))+(28*x(12*(j-1)+2))+(31*x(12*(j-1)+3)))/90 ;
        end
        return
    else

    if (season == 'FMA')
        for j=1:length(x)/12
            seasonal_mean(j) = ((28*x(12*(j-1)+2))+(31*x(12*(j-1)+3))+(30*x(12*(j-1)+4)))/89 ;
        end
        return
    else

    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end


elseif (size(season) == [1 2])

    if (season == 'DJ')
        for j=2:length(x)/12-1
            seasonal_mean(j) = ((31*x(12*(j-1)+12))+(31*x(12*(j)+1)))/62 ;
        end
        seasonal_mean(1) = seasonal_mean(2);
        return
    end
    if (season == 'JJ')
        for j=1:length(x)/12
            seasonal_mean(j) = ((30*x(12*(j-1)+6))+(31*x(12*(j-1)+7)))/61 ;
        end
        return
    end
end
