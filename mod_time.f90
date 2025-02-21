module mod_time

  use mod_param
  implicit none
  integer :: di33
  integer , dimension(13) :: dfc
  character(9) , dimension(12) :: emonths , months
  character(11) , dimension(7) :: eweekd , weekd
  data months/'Gennaio' , 'Febbraio' , 'Marzo' , 'Aprile' , 'Maggio' ,   &
              'Giugno' , 'Luglio' , 'Agosto' , 'Settembre' , 'Ottobre' , &
              'Novembre' , 'Dicembre'/
  data emonths/'January' , 'February' , 'March' , 'April' , 'May' ,   &
               'June' , 'July' , 'August' , 'September' , 'October' , &
               'November' , 'December'/
  data eweekd/'Monday' , 'Tuesday' , 'Wednesday' , 'Thursday' , 'Friday' ,&
              'Saturday' , 'Sunday'/
  integer , dimension(100) :: iflg
  data iflg(1)/0/            ! log level
  data iflg(2)/700/          ! X size grafici
  data iflg(3)/1/            ! boxplot should or should not draw the frame
  data iflg(4)/0/            ! Plot 1D
  data iflg(5)/0/            ! Dominio 0-3
  data iflg(6)/1/            ! Mese Corrente
  data iflg(7)/1/            ! Check Landuse
  data iflg(8)/0/            ! Logo sui plot
  data iflg(9)/0/            ! Stile HTML
  data iflg(10)/2/           ! Stile con cui produce le date
  data iflg(11)/-1/          ! confini politici (no longer used)
  data iflg(12)/1/           ! colore delle scritte dei grafici
  data iflg(13)/1/           ! colore degli assi nei plot 1-d
  data iflg(14)/4/           ! colore delle freccette del vento
  data iflg(15)/0/           ! Se > 0 nella sub. meteogrammi setta l'asse x
  data iflg(16)/0/           ! colore dei confini politici
  data iflg(17)/-1/          ! se >= 0 disegna i confini di default gks ncar
  data iflg(18)/1/           ! unita' di misura u-v (1=m/s 2=Km/h)
  data iflg(19)/1/           ! colore delle isolinee (routine mm5contorni)
  data iflg(20)/0/           ! Numero di curve di livello
  data iflg(21)/2/           ! Dimensioni delle scritte (1,2,3 della wtstr)
  data iflg(22)/1/           ! coord.grafiche (1=x,y 2=lat-lon 3=correnti)
  data iflg(23)/0/           ! primo indice delle matrici (0=lon altro=lat)
  data iflg(24)/0/           ! step tra i livelli degli automi cellulari
  data iflg(25)/1/           ! colore sfondo dei grafici: 0=nero altro=bianco
  data iflg(26)/18/          ! numero di colori
  data iflg(27)/1/           ! palette di colori
  data iflg(28)/0/           ! contorno ai pallocchi
  data iflg(29)/1/           ! La label bar di mm5colormap
  data iflg(30)/0/           ! se = 0 viacolvento cancella i singoli frames
  data iflg(31)/0/           ! label asse x ( 1 = ore 2 = mesi 3 = nulla )
  data iflg(32)/0/           ! Numero di colori nei plot 2d
  data iflg(33)/1/           ! Ogni quanti punti griglia fa le freccette
  data iflg(34)/875949887/   ! Random seed for gauss and flat distribution
  data iflg(35)/0/           ! Qnorm - number of bin to be used
  data iflg(36)/0/           ! label asse y ( 3 = nulla )
  data iflg(37)/1/           ! Primo colore per disegnare le palette
  data iflg(38)/1/           ! Ogni quanti colori per disegnare le palette
  data iflg(39)/1/           ! Scrive il massimo e il minimo - subr. contorni
  data iflg(40)/-1/          ! Colore Plotta
  data iflg(41)/-13/         ! Tipo di carattere vedi man gstxfp
  data iflg(42)/20/          ! Colore delle freccette (subr. freccetelle)
  data iflg(43)/0/           ! Colore sfondo dei plot 1-d
  data iflg(44)/0/           ! Used in plotta package - internal only
  data iflg(45)/1/           ! boundary routine also write main town if > 0
  data iflg(46)/3/           ! Visulaizzaplot format (0=bmp, 1=gif, 2=tiff)
  data iflg(47)/620/         ! Size delle finestre negli script html
  data iflg(48)/0/           ! nei plot mm5 se = 1 disegna l'italia
                             !              se = 2 disegna le regioni
                             !              se = 3 disegna le provincie
  data iflg(49)/0/           ! Se = 1 disegna regioni Se = 2 provincie
  data iflg(50)/0/           ! Internal use
  data iflg(51)/0/           ! Se = 1 label orizzont  Se = 2 label vertic
  data iflg(52)/0/           ! Stile delle freccette (routine boxarrow)
  data iflg(53)/-1/          ! Histogram color (routine mvbookplot)
  data iflg(54)/20/          ! Numero di colori della palette corrente
  data iflg(55)/0/           ! Internal use only - plottaframe, plottavw
  data iflg(56)/-1/          ! Size delle label - plotpalette, boxlegend
  data iflg(57)/-1/          ! secolo - subroutine datafrommm5index ecc.
  data iflg(58)/1/           ! controllo video. <>0 = ON
  data iflg(59)/0/           ! Visulaizzaplot behaviour 0=vis ; 1=salva
  data iflg(60)/1/           ! Used by somepoch to define neural initializ.
  data iflg(61)/1/           ! Used by somepoch to define neural distance
  data iflg(62)/0/           ! Used by openmuseofile
  data iflg(63)/0/           ! Used by openmuseofile (internal)
  data iflg(64)/0/           ! Used by boxplot
  data iflg(65)/0/           ! Used by chymdownscaling
  data iflg(66)/0/           ! Background color of html table
  data iflg(67)/3/           ! Histogram frame (subroutine mvbookplot)
  data iflg(68)/1/           ! Histogram statistic (subroutine mvbookplot)
  data iflg(69)/0/           ! Plotta fill
  data iflg(70)/6/           ! Logical unit used by CHyM Model for log
  data iflg(71)/0/           ! CHyM Standard domain
  data iflg(72)/0/           ! Calendar

  contains
  subroutine gmafromindex(dindex,ih,id,im,iy)
    implicit none
    integer :: id , ih , im , dindex , iy
    intent (in) dindex
    intent (out) id , ih , im , iy
    character(10) :: tmps
    write (tmps,'(i10)',err=100) dindex
    read (tmps,'(i4,3i2)',err=100) iy , im , id , ih
    return
 100  iy = -9999
      im = -9999
      id = -9999
      ih = -9999
  end subroutine gmafromindex
  integer function increasetime(dindex)
    implicit none
    integer :: dindex
    intent (in) dindex
    integer :: id , ih , iii , ilm , im , iy
    character(10) :: tmps
    if ( dindex==1999123123 .or. (trim(calendario) =='360_day' .and. dindex==1999123023) ) then
      increasetime = 2000010100
      return
    end if
    write (tmps,'(i10)',err=100) dindex
    read (tmps,'(i4,3i2)',err=100) iy , im , id , ih
    if ( ical==1 ) then
      ilm = monthlen(im,2001)
    else if ( ical==2 ) then
      ilm = 30
    else
      ilm = monthlen(im,iy)
    end if
    ih = ih + 1
    if ( ih==24 ) then
      ih = 0
      id = id + 1
      if ( id>ilm ) then
        id = 1
        im = im + 1
        if ( im==13 ) then
          im = 1
          iy = iy + 1
        end if
      end if
    end if
    if ( ih<0 .or. ih>23 .or. im<1 .or. im>12 .or. id<1 .or. id>ilm .or.    &
         iy<1900 .or. iy>2099 ) then
      increasetime = -9999
    else
      increasetime = ih + id*100 + im*10000 + iy*1000000
    end if
    return
 100  increasetime = -9999
  end function increasetime
  subroutine dataorafromday(ih,id,im,iy,cdata)
    implicit none
    character(len=*) :: cdata
    integer :: id , ih , im , iy
    intent (in) id , ih , im , iy
    intent (out) cdata
    call datafromindex(indexofyear(id,im,iy),iy,cdata)
    if ( ih>=10 ) then
      write (cdata,'(a,i2)') trim(cdata)//' h: ' , ih
    else
      write (cdata,'(a,i1)') trim(cdata)//' h: 0' , ih
    end if
    call no2space(cdata)
    if ( iflg(10)==4 ) cdata = trim(cdata)//' UGT'
    if ( iflg(10)==5 ) cdata = trim(cdata)//' CEST'
  end subroutine dataorafromday
  subroutine datafromindex(dindex,anno,cdata)
    implicit none
    integer :: anno , dindex
    character(len=*) :: cdata
    integer :: id , im , iy
    iy = i4digityear(anno)
    id = -1
    im = -1
    call dayfromindex(dindex,id,im,iy)
    if ( id>0 .and. im>0 ) then
      if ( iflg(10)==1 ) then
        write (cdata,'(i2,1x,a9,i5)') id , months(im) , iy
      else if ( iflg(10)==2 ) then
        write (cdata,'(a,i2,a,i5)') trim(emonths(im))//' ' , id , ', ' , iy
      else if ( iflg(10)>=3 .and. iflg(10)<=5 ) then
        write (cdata,'(a11,1x,a9,1x,i2,i5)') eweekd(iweekday(id,im,iy)) , &
                                             emonths(im) , id , iy
      else if ( iflg(10)==6 ) then
        write (cdata,'(a,i2,a,i5)') emonths(im)(1:3)//' ' , id , ', ' , iy
      else if ( iflg(10)==7 ) then
        write (cdata,'(a4,i2)') emonths(im)(1:3)//' ' , id
      else if ( iflg(10)==8 ) then
        write (cdata,'(a3,i2)') emonths(im)(1:3) , id
      else
        write (cdata,'(a11,1x,i2,1x,a9,i5)') weekd(iweekday(id,im,iy)) , id , &
                                             months(im) , iy
      end if
      call no2space(cdata)
      call noinspace(cdata)
    end if
  end subroutine datafromindex
  integer function monthlen(im,iy)
    implicit none
    integer :: im , iy
    intent (in) im , iy
    if ( di33/=33 .or. iy/=oldyear ) then
      oldyear = iy
      call definizionedimvtime
    end if
    monthlen = mesi(im)
  end function monthlen
  subroutine dayfromindex(dindex,id,im,iy)
    implicit none
    integer :: id , im , dindex , iy
    intent (in) dindex , iy
    intent (out) id , im
    integer :: i
    if ( di33/=33 .or. iy/=oldyear ) then
      oldyear = iy
      call definizionedimvtime
    end if
    do i = 1 , 12
      if ( dindex>dfc(i) .and. dindex<=dfc(i+1) ) then
        im = i
        id = dindex - dfc(i)
        return
      end if
    end do
    write (6,'(a,i10)') ' Bad index passed to dayfromindex ' , dindex
  end subroutine dayfromindex

  subroutine noinspace(a)
    implicit none
    character(len=*) , intent(inout) :: a
    a = adjustl(a)
  end subroutine noinspace
  subroutine no2space(a)
    implicit none
    character(len=*) :: a
    intent (inout) a
    integer :: i
    no2space_loop: &
    do
      do i = 1 , len_trim(a) - 1
        if ( ichar(a(i:i))==32 .and. ichar(a(i+1:i+1))==32 ) then
          a = a(1:i)//a(i+2:)
          cycle no2space_loop
        end if
      end do
      exit
    end do no2space_loop
  end subroutine no2space
  integer function indexofyear(id,im,iy)
    implicit none
    integer :: id , im , iy
    intent(in) :: id , im , iy
    indexofyear = index1d(id,im,iy)
  end function indexofyear
  integer function index1d(id,im,iy)
    implicit none
    integer :: id , im , iy
    intent (in) id , im , iy
    if ( di33/=33 .or. iy/=oldyear ) then
      oldyear = iy
      call definizionedimvtime
    end if
    index1d = dfc(im) + id
  end function index1d
  integer function iweekday(id,im,iy)
    implicit none
    integer :: id , im , iy
    intent (in) id , im , iy
    integer :: igiorni
    if ( di33/=33 .or. iy/=oldyear ) then
      oldyear = iy
      call definizionedimvtime
    end if
    if ( mod(iy,4)==0 ) then
      igiorni = id + dfc(im) + iy*365 + iy/4 + 3
    else
      igiorni = id + dfc(im) + iy*365 + iy/4 + 4
    end if
    iweekday = mod(igiorni,7) + 1
  end function iweekday
  subroutine definizionedimvtime ! internal mvlib use
    implicit none
    integer :: i
    if ( mod(oldyear,4)==0 ) then
      mesi(2) = 29
    else
      mesi(2) = 28
    end if
    dfc(1) = 0
    do i = 2 , 13
      dfc(i) = dfc(i-1) + mesi(i-1)
    end do
    di33 = 33
  end subroutine definizionedimvtime
  integer function i4digityear(anno)
    implicit none
    integer :: anno
    intent (in) anno
    i4digityear = anno
    if ( i4digityear<1000 ) then
      if ( iflg(57)>0 ) then
        i4digityear = i4digityear + iflg(57)
      else if ( i4digityear<=50 ) then
        i4digityear = i4digityear + 2000
      else
        i4digityear = i4digityear + 1900
      end if
    end if
  end function i4digityear

end module mod_time
