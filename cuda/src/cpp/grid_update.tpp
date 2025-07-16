template<typename prcs_type>
void update ()
{
    for ( int i_soln = 0; i_soln < this->N_soln; i_soln++ )
    {
        run_time_vector<prcs_type> & S = this->Vec_soln.at(i_soln);
        run_time_vector<prcs_type> & R = this->Vec_rths.at(i_soln);

        for ( int i = 0; i < G_size_x * G_size_y; i++ )
        {
            S(i) = S(i) + R(i);
            R(i) = 0.;
        }
    }
} // update ()


template<typename prcs_type>
void update_3C ()
{
    for ( int i_soln = 0; i_soln < this->N_soln; i_soln++ )
    {
        run_time_vector<prcs_type> & S = this->Vec_soln.at(i_soln);
        run_time_vector<prcs_type> & R = this->Vec_rths.at(i_soln);
        run_time_vector<prcs_type> & C = this->Vec_cpst.at(i_soln);

        for ( int i = 0; i < G_size_x * G_size_y; i++ )
        {
            prcs_type Y = R(i) + C(i);
            prcs_type T = S(i) + Y;
            C(i) = T - S(i);
            C(i) = Y - C(i);
            S(i) = T;

            R(i) = 0.;
        }
    }
} // update_3C ()


template<typename prcs_type>
void update_3R ()
{
    for ( int i_soln = 0; i_soln < this->N_soln; i_soln++ )
    {
        run_time_vector<prcs_type> & S = this->Vec_soln.at(i_soln);
        run_time_vector<prcs_type> & R = this->Vec_rths.at(i_soln);

        for ( int i = 0; i < G_size_x * G_size_y; i++ )
        {
            prcs_type Y = R(i);
            prcs_type T = S(i) + Y;
            R(i) = T - S(i);
            R(i) = Y - R(i);
            S(i) = T;
        }
    }
} // update_3R ()


template<typename prcs_type>
void update_3C_fast2sum ()
{
    for ( int i_soln = 0; i_soln < this->N_soln; i_soln++ )
    {
        run_time_vector<prcs_type> & S = this->Vec_soln.at(i_soln);
        run_time_vector<prcs_type> & R = this->Vec_rths.at(i_soln);
        run_time_vector<prcs_type> & C = this->Vec_cpst.at(i_soln);

        for ( int i = 0; i < G_size_x * G_size_y; i++ )
        {
            prcs_type Y = R(i) + C(i);
            fast2sum<prcs_type> ( S(i) , Y , S(i) , C(i) );
            R(i) = 0.;
        }
    }
} // update_3C_fast2sum ()


template<typename prcs_type>
void update_3R_fast2sum ()
{
    for ( int i_soln = 0; i_soln < this->N_soln; i_soln++ )
    {
        run_time_vector<prcs_type> & S = this->Vec_soln.at(i_soln);
        run_time_vector<prcs_type> & R = this->Vec_rths.at(i_soln);

        for ( int i = 0; i < G_size_x * G_size_y; i++ )
        {
            fast2sum<prcs_type> ( S(i) , R(i) , S(i) , R(i) );
        }
    }
} // update_3R_fast2sum ()


template<typename prcs_type>
void update_6C_slow2sum ()
{
    for ( int i_soln = 0; i_soln < this->N_soln; i_soln++ )
    {
        run_time_vector<prcs_type> & S = this->Vec_soln.at(i_soln);
        run_time_vector<prcs_type> & R = this->Vec_rths.at(i_soln);
        run_time_vector<prcs_type> & C = this->Vec_cpst.at(i_soln);

        for ( int i = 0; i < G_size_x * G_size_y; i++ )
        {
            prcs_type Y = R(i) + C(i);
            slow2sum<prcs_type> ( S(i) , Y , S(i) , C(i) );
            R(i) = 0.;
        }
    }
} // update_6C_slow2sum ()


template<typename prcs_type>
void update_6R_slow2sum ()
{
    for ( int i_soln = 0; i_soln < this->N_soln; i_soln++ )
    {
        run_time_vector<prcs_type> & S = this->Vec_soln.at(i_soln);
        run_time_vector<prcs_type> & R = this->Vec_rths.at(i_soln);

        for ( int i = 0; i < G_size_x * G_size_y; i++ )
        {
            slow2sum<prcs_type> ( S(i) , R(i) , S(i) , R(i) );
        }
    }
} // update_6R_slow2sum ()