    PROGRAM SIMUL_IPH2_HORA

    implicit none


    !===============================================================================
    ! DECLARAÇÃO DE VARIÁVEIS
    !===============================================================================
    integer, parameter :: ndim = 9           !Número de parâmetros do modelo hidrológico
    integer :: i, j, k, eof, ND
    integer, dimension(100000) :: ano, mes, dia, hora
    real*8, dimension(100000) :: ETp, CMB, Qmont, Qexut, peso, Qmod
    real*8 :: Ainc, Atot, estats(9), P(ndim)
    character :: arqresult*100, arqecv*100, frt*7


    !===============================================================================
    ! INICIALIZAÇÕES
    !===============================================================================
    frt = "9f12.6"




    !===============================================================================
    ! LEITURA DE ENTRADAS
    !===============================================================================
    !Opções da rodada
    open(unit=11, file='config_iph2.txt', action='read')
    read(11,*) (P(i), i = 1, ndim)
    read(11,'(a)') arqecv
    read(11,'(a)') arqresult
    read(11,*) Atot, Ainc
    close(11)
    
    write(*,'(a,2(f9.3, a))') "Area Total = ", Atot, " km2;    Area Incremental = ", Ainc, " km2."
    write(*,'(a,20f12.6)') "Parametros: ", (P(i), i = 1, ndim)


    !Contando linhas (quantidade de dados) do arquivo 'ecv.txt'
    open(unit=10, file=arqecv, action='read')
    do i = 1, 100000
        read(10,*,iostat=eof) ano(i), mes(i), dia(i), hora(i), ETp(i), CMB(i), Qmont(i), Qexut(i), peso(i)
        if (eof < 0) exit
    end do
    close(10)
    ND = i - 1
    write(*,'(/,a,i6.6,a)') "Armazenou ", ND, " registros horários."




    !===============================================================================
    ! GRAVANDO SAÍDAS
    !===============================================================================

    !Gravando resultados da calibração
    open(unit = 12, file = adjustl(arqresult), action = 'write')
    call FObj(ETp, CMB, Qmont, Qexut, peso, ND, Ainc, P, estats, Qmod)
    write(12,'(a1, 6f9.3, 3f10.6)')'#', (estats(i), i = 1, 9)
    
    !Vazão modelada
    do i = 1, ND
        write(12,'(i4,3i3.2,f8.2)') ano(i), mes(i), dia(i), hora(i), Qmod(i)
    end do
    
    write(*,*)'Simulacao HORARIA concluida!'
    
    
    END PROGRAM SIMUL_IPH2_HORA




    !===============================================================================
    ! Estatísticas: Parâmetros estatísticos da qualidade da vazão simulada
    !===============================================================================
    SUBROUTINE FObj(ETp, CMB, Qmont, Qexut, peso, N, Area, Params, func, Qmod)
    implicit none
    
    integer :: i, N
    real*8, dimension(100000) :: ETp, CMB, Qmont, Qexut, Qmod, peso
    real*8 :: Params(9), Area, func(9), media(3), desvp(3), sp, res
    
    !Gerando série de vazão modelada
    call IPH2hor(ETp, CMB, Qmont, Qexut, N, Area, Params, Qmod)
    
    !Calculando valor média da vazão na exutória, vazão modelada e dos resíduos
    media = 0.d0
    do i = 1, N
        media(1) = media(1) + Qexut(i) * peso(i) !acumulando valores de vazão observada
        media(2) = media(2) + Qmod(i) * peso(i)  !acumulando valores de vazão modelada
        media(3) = media(3) + (Qmod(i) - Qexut(i)) * peso(i)  !acumulando valores do resíduo
    end do
    sp = sum(peso)
    media(1:3) = media(1:3)/sp !média ponderada
    
    !Calculando desvio-padrão das séries de vazão e de resíduo
    !Faz o somatório das séries de erro e variância
    desvp = 0.d0
    func  = 0.d0
    do i = 1, N
        res = Qmod(i) - Qexut(i) !Resíduo sem a ponderação
        desvp(1) = desvp(1) + ((Qexut(i) - media(1))**2) * peso(i)
        desvp(2) = desvp(2) + ((Qmod(i) - media(2))**2) * peso(i)
        desvp(3) = desvp(3) + ((res - media(3))**2) * peso(i)
        func(1) = func(1) + dabs(res) * peso(i)
        func(2) = func(2) + ((Qmod(i) - Qexut(i))**2) * peso(i)
        if (res > 0) then
            func(3) = func(3) + res * peso(i)
            func(4) = func(4) + peso(i)
        end if
        if (res < 0) then
            func(5) = func(5) + res * peso(i)
            func(6) = func(6) + peso(i)
        end if
        func(7) = func(7) + (Qexut(i) - media(1)) * (Qmod(i) - media(2)) * peso(i)
        func(8) = func(8) + ((Qexut(i) - media(1))**2) * peso(i)
    end do
    desvp(1:3) = dsqrt(desvp(1:3)/sp)
    
    !Contabilizando estatísticas
    func(1) = func(1)/sp                           !Erro Absoluto Médio (m3/s)
    func(8) = 1.d0 - func(2)/func(8)               !Coeficiente de Nash-Sutcliffe (adim.)
    func(2) = dsqrt(func(2)/sp)                    !Raiz do Erro Quadrático Médio (m3/s)
    func(3) = func(3)/func(4)                      !Erro Positivo Médio (m3/s)
    func(4) = func(4)*100.d0/sp                    !Frequência de Erros Positivos (%)
    func(5) = func(5)/func(6)                      !Erro Negativo Médio (m3/s)
    func(6) = func(6)*100.d0/sp                    !Frequência de Erros Negativos (%)
    func(7) = func(7)/(desvp(1)*desvp(2)*sp)       !Correlação linear de Pearson (adim.)
    func(9) = 0.3d0 * dabs(media(3))/media(1) + 0.7d0 * desvp(3)/desvp(1)    !Coeficiente sem nome (adim.)
    
    END SUBROUTINE FObj




    !===============================================================================
    ! Modelo: Instituto de Pesquisas Hidráulicas II - IPH 2
    !===============================================================================
    SUBROUTINE IPH2hor(ETp, CMB, Qin, Qout, N, Area, Params, Qbac)
    implicit none
    
    integer :: i, j, N
    real*8, dimension(100000) :: ETp, CMB, Qin, Qout, Qbac
    real*8 :: Params(9), Area, conv
    
    integer :: KT, NH, NDT
    real*8 :: RMAX, Io, Ib, H, ALF, Ks, Kb, Ainp, Kprop
    real*8 :: P, E, EP, R, ER, S, Smax, RI, T, RD, BI, AI, AT, BT
    real*8 :: AIL, BIL, ATL, BTL    
    real*8 :: AT1, Par, CR, S1, RI1, VE, VI, VP, SX, ATX, RAUX, VAUX
    real*8 :: PV(100), HIST(100), QS, QT, TCI
    
    !Repassando parâmetros do modelo para variáveis
    RMAX = Params(1)
    Io   = Params(2)
    Ib   = Params(3) * Io
    H    = Params(4)
    ALF  = Params(5)
    Ks   = dexp(-1.d0/Params(6))
    Kb   = dexp(-1.d0/Params(7))
    Ainp = Params(8)
    NH   = nint(Params(9))
    
    !Constantes do modelo computadas sobre os valores dos parâmetros
    Smax = -Io / dlog(H)
    BI   = Io / dlog(H) / (Io - Ib)
    AI   = -Io * BI
    AT   = 0.d0
    BT   = -Io / Ib / dlog(H)
    AIL  = -AI / BI
    BIL  = 1.d0 / BI
    ATL  = 0.d0
    BTL  = 1.d0 / BT
    conv = Area / 3.6d0
    
    !Inicializações
    S  = 0.5d0 * Smax
    R  = 0.5d0 * RMAX !0.d0
    RI = AIL + BIL * S
    PV = 0.d0
    QT = max((Qout(1) - Qin(1)), 0.d0) / conv
    QS = 0.d0
    HIST = 0.d0
    HIST(1:NH) = 1.d0/dfloat(NH) !Bacia retangular pro histograma tempo-área

    !Rodando modelo IPH2 para a série de dados
    do i = 1, N
    
        P = CMB(i)
        E = ETp(i)

    ! 1   - Perdas por evaporação e interceptação
    !------------------------------------------------------------------------------------------------------------------+
    !   A evaporação potencial é retirada da precipitação quando for inferior a esta, e em caso contrário, a evaporação|
    !potencial não satisfeita é atendida pelo reservatório de interceptação (cobertura vegetal e depressões). Quando   |
    !este último reservatório está totalmente esgotado, o déficit de evaporação potencial passa a ser atendido pela    |
    !água contida no solo, através da relação linear entre a evapotranspiração e a porcentagem de umidade do solo, dado|
    !por:                                                                                                              |
    !E(t) = EP(t) * ( S(t) / Smax ),    onde E(t) é a evapotranspiração da superfície no tempo t; EP(t) é a evapotrans-|
    !piração potencial; S(t) é o estado de umidade da camada superior do solo.                                         |
    !   Quando a precipitação é maior que a evaporação potencial, a diferença é retida por interceptação até que sua   |
    !capacidade máxima Rmax seja satisfeita.                                                                           |
    !------------------------------------------------------------------------------------------------------------------+
        if (P < E) then
            EP = E - P
            P  = 0.0

            if (EP <= R) then
                R = R - EP

            else
                EP = EP - R
                R  = 0.0
                ER = EP * S/Smax
                S  = S - ER
                if (S < 0.0) then
                    ER = ER + S
                    S  = 0.0
                end if

                RI = AIL + BIL*S
                T  = ATL + BTL*S
                ER = ER + P

            end if

        else
            P  = P - E
            ER = E
            RD = RMAX - R

            if (P <= RD) then
                R = R + P
                P = 0.0
            else
                P = P - RD
                R = RMAX
            end if

        end if

    ! 2   - Separação dos volumes
    !------------------------------------------------------------------------------------------------------------------+
    !   A parcela de precipitação resultante pode gerar escoamento superficial ou infiltrar no solo, entretanto a      |
    !parcela de água que precipita sobre áreas impermeáveis da bacia gera escoamento superficial direto sem que ocorra |
    !infiltrações. A fração da bacia coberta por áreas impermeáveis é dada pelo parâmetro AIMP.                        |
    !   Da parcela de água que precipitou sobre a área permeável da bacia é necessário calcular o volume percolado para|
    !o aquífero e o volume que gera escoamento superficial. Pela equação da continuidade tem-se o seguinte:            |
    !dS/dt = I(t) - T(t),    sendo S o estado de umidade do solo (mm), I(t) a infiltração e T(t) a percolação.         |
    !   A infiltração é contabilizada pela equação de Horton e a percolação por uma fórmula proposta por Berthelot:    |
    !I(t) = Ib + (Io - Ib)*(h**dt),    T(t) = Ib*(1 - h**dt),    onde Ib é a capacidade de infiltração quando o solo   |
    !está saturado, Io é a capacidade de infiltração do solo quando a umidade é So, h = e**(-k), e k é um parâmetro que|
    !caracteriza o decaimento da curva exponencial de infiltração e depende das caracteristicas do solo.               |
    !------------------------------------------------------------------------------------------------------------------+
        AT1 = 1.d0
        Par = P

        if (P < RI) then
            CR  = (P/RI)**2 / ((P/RI)+ALF)
            P   = P*(1.d0 - CR)
            S1  = (S*(2.d0 - 1.d0/BT) + 2.d0*P) / (2.d0 + 1.d0/BT)
            RI1 = AIL + BIL*S1

            if (P < RI1) then
                T  = ATL + BTL*S1
                VE = 0.d0
                VI = P
            else
                SX   = AI + BI*P
                ATX  = 2.d0*BT*(SX-S) / (2.d0*P*BT + 2.d0*AT - SX - S) !AT = 0
                AT1  = AT1 - ATX
                RAUX = P
                VAUX = P*ATX

                RI1 = Ib + (RAUX-Ib) * H**AT1
                S1  = AI + BI*RI1
                T   = ATL + BTL*S1
                VI  = Ib*AT1 + (RAUX-Ib)*(H**AT1 -1.d0)/dlog(H) + VAUX
                VE  = P*AT1 - VI + VAUX
            end if

        else
            RAUX = RI
            VAUX = 0.d0

            RI1 = Ib + (RAUX-Ib) * H**AT1
            S1  = AI + BI*RI1
            T   = ATL + BTL*S1
            VI  = Ib*AT1 + (RAUX-Ib)*(H**AT1 -1.)/dlog(H) + VAUX
            VE  = P*AT1 - VI + VAUX
        end if

        VP = S - S1 + VI
        VE = VE*(1.d0-Ainp) + Par*Ainp

    ! 3   - Propagação dos escoamentos
    !------------------------------------------------------------------------------------------------------------------+
    !   O escoamento superficial é propagado pelo modelo Clark o qual utiliza um histograma tempo-área para simular o  |
    !deslocamento da água ao longo da bacia e o método de reservatório linear para o efeito de atenuação. Para o esco- |
    !amento subterrâneo apenas o método de reservatório linear é utilizado. A matriz do histograma, HTA, é na realidade|
    !dois vetores. Na primeira coluna devem estar os pesos do histograma para cada seção e na segunda coluna os volumes|
    !do escoamento superficial acumulados para cada seção da bacia.                                                    |
    !------------------------------------------------------------------------------------------------------------------+
        do KT = 1, NH
            PV(KT) = PV(KT) + VE*HIST(KT)
        end do
        VE = PV(1)

        do KT = 1, NH-1
            PV(KT) = PV(KT+1)
        end do
        PV(NH) = 0.d0

        QS = QS*Ks + VE*(1.d0-Ks)
        QT = QT*Kb + VP*(1.d0-Kb)
        S  = S1
        RI = RI1
        
        Qbac(i) = (QS + QT) * conv
    
    end do

    END SUBROUTINE IPH2hor
