void NewtonRaphson(std::vector<Node*> &nodes, std::vector<Element*> &elements, std::vector<Physical*> &physicals, std::vector<Parameter*> &parameters,
					std::map<Node*, std::vector<double> > &solution, std::vector<double> &theta_k,gmm::row_matrix< gmm::wsvector<double> > Tmp,
					std::vector<double> &qext,FemFlag method,map<int, double> linesRegion, map<Node*, Node*> &NodesCorresp,
					std::vector<double> &delta_theta_k)
{

	std::vector<double> qint(nodes.size());
    gmm::csr_matrix<double> KT;

	Internal_flux(std::vector<double> &theta_k, std::map<int,std::vector<double> > &surfaceRegion, std::vector<Element*> &elements, std::vector<double> &qint);

    Tangent_Stiffness_Matrix(std::vector<double> &theta_k,std::map<int,std::vector<double> > &surfaceRegion, std::vector<Element*> &elements,
                             gmm::csr_matrix<double> &KT,map<Node*, Node*> &NodesCorresp);

    //Conditions aux limites sur KT
    if(method == DIRICHLETFLAG)
    {
        for(unsigned int i = 0; i < elements.size(); i++)
        {
            if(elements[i]->type == 1)//If line
            {
                if(linesRegion.count(elements[i]->region) == 1)//If linesRegion contains elements[i]->region
                {
                    for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
                    {
                        for(unsigned int k = 0; k < nodes.size(); k++)
                        {
                            if(k == elements[i]->nodes[j]->num-1)
                            {
                                KT(elements[i]->nodes[j]->num-1, k) = 1;
                            }
                            else
                            {
                                KT(elements[i]->nodes[j]->num-1, k) = 0;
                            }
                        }
                    }
                }
            }
        }

        //Dirichlet sur f

        for(unsigned int i = 0; i < elements.size(); i++)
        {
            if(elements[i]->type == 1)//If line
            {
                if(linesRegion.count(elements[i]->region) == 1)//If linesRegion contains elements[i]->region
                {
                    for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
                    {
                        f[elements[i]->nodes[j]->num-1] = 0;
                    }
                }
            }
        }
    }

    std::vector<double> RHS(nodes.size()); //RHS of Newton Raphson system
    for (i=0 ; i<nodes.size(); i++)
        RHS(i) = qint(i)-qext(i);

    //System
    gmm::lu_solve(KT, delta_theta_k, RHS);

    for (i=0; i<nodes.size(); i++)
        theta_k(i)+=delta_theta_k(i);


}

void Internal_flux(std::vector<double> &theta_k, std::map<int,std::vector<double> > &surfaceRegion, std::vector<Element*> &elements, std::vector<double> &qint)
{

    gmm::dense_matrix<double> alpha(2,2); // matrice alpha

    alpha(0,0) = surfaceRegion[elements[i]->region][1];
    alpha(0,1) = surfaceRegion[elements[i]->region][2];
    alpha(1,0) = surfaceRegion[elements[i]->region][3];
    alpha(1,1) = surfaceRegion[elements[i]->region][4];

    gmm::dense_matrix<double> beta(2,2); // matrice beta pour la conductivité

    beta(0,0) = surfaceRegion[elements[i]->region][5];
    beta(0,1) = surfaceRegion[elements[i]->region][6];
    beta(1,0) = surfaceRegion[elements[i]->region][7];
    beta(1,1) = surfaceRegion[elements[i]->region][8];

	for(unsigned int i = 0; i < elements.size(); i++)//loop over the elements
    {
        if(elements[i]->type == 2)//If triangle
        {
            Node *n1 = elements[i]->nodes[0];
            Node *n2 = elements[i]->nodes[1];
            Node *n3 = elements[i]->nodes[2];
            double x1 = n1->x;
            double y1 = n1->y;
            double x2 = n2->x;
            double y2 = n2->y;
            double x3 = n3->x;
            double y3 = n3->y;
            double qint_I; //Contribution of a given node to qint

            std::vector<double> theta_nodes(3);//Temperature at nodes
            theta_nodes(0) = theta_k[n1->num-1];
            theta_nodes(1) = theta_k[n2->num-1];
            theta_nodes(2) = theta_k[n3->num-1];
            double thetam = (theta_nodes(0)+theta_nodes(1)+theta_nodes(2))/3;//Average temperature on the element
            double detJ = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);


            gmm::dense_matrix<double> kappa(2,2); //Conductivity tensor at given temperature

            for (unsigned int i = 0; i<2;i++)
            {
            	for(unsigned int j = 0 ; j<2; j++)
            	{
            		kappa(i,j) = alpha (i,j)*thetam+beta(i,j);
            	}
            }


            gmm::dense_matrix<double> inv_J(2,2); //Inverse of the Jacobian Matrix
            inv_J(0,0) = 1/detJ*(y3-y1);
            inv_J(0,1) = 1/detJ*(y1-y2);
            inv_J(1,0) = 1/detJ*(x3-x1);
            inv_J(1,1) = 1/detJ*(x2-x1);

            gmm::dense_matrix<double> grad_phi(2,3);
            grad_phi(0,0) = -1;
            grad_phi(1,0) = -1;
            grad_phi(0,1) = 1;
            grad_phi(1,1) = 0;
            grad_phi(0,2) = 0;
            grad_phi(1,2) = 1;

            //Computation of qint

            std::vector<double> tmp_vec(2);//Temporary vector in the computation of qint
            std::vector<double> tmp_vec_2(2);//Temporary vector in the computation of qint

            tmp_vec(0)=0;
            tmp_vec(1)=0;

            for (unsigned int i = 0;i<2;i++)
            {
            	for(unsigned int j = 0; j<3 ; j++)
            	{
            		tmp_vec(0) = tmp_vec(0)+theta_nodes(j)*grad_phi(0,j);
            		tmp_vec(1) = tmp_vec(1)+theta_nodes(j)*grad_phi(1,j);
            	}
            	for(unsigned int j = 0; j<2; j++){
            		for(unsigned int k = 0; k<2 ; k++)
            		{
            			tmp_vec_2(j) = inv_J(j,k)*tmp_vec(k);
            		}
            	}
            	for(unsigned int j = 0; j<2; j++){
            		for(unsigned int k = 0; k<2 ; k++)
            		{
            			tmp_vec(j) = kappa(j,k)*tmp_vec_2(k);
            		}
            	}
            }

            for(unsigned int i = 1; i<4 ; i++)//For each element, there are three contributions to qint at distinct nodes
            {
            	for(unsigned int j = 0; j<2; j++)
                {
            		for(unsigned int k = 0; k<2 ; k++)
            		{
            			tmp_vec_2(j) = inv_J(j,k)*grad_phi(k,i);
            		}
            	}
            }

            	qint_I = -detJ/2*(tmp_vec(0)*tmp_vec_2(0)+tmp_vec(1)*tmp_vec_2(1));

            	if(i == 1)
                    qint[n1->num-1] += qint_I;
                if(i == 2)
                    qint[n2->num-1] += qint_I;
                if(i == 3)
                    qint[n3->num-1] += qint_I;
        }


    }


}


void Tangent_Stiffness_Matrix(std::vector<double> &theta_k,std::map<int,std::vector<double> > &surfaceRegion, std::vector<Element*> &elements, gmm::csr_matrix<double> &KT,map<Node*, Node*> &NodesCorresp)
{
    //K matrix
    gmm::row_matrix< gmm::wsvector<double> > Tmp(nodes.size(), nodes.size());

    for(unsigned int i = 0; i < elements.size(); i++)
    {
        if(elements[i]->type == 2)//If triangle
        {
            Node *n1 = elements[i]->nodes[0];
            Node *n2 = elements[i]->nodes[1];
            Node *n3 = elements[i]->nodes[2];
            double x1 = n1->x;
            double y1 = n1->y;
            double x2 = n2->x;
            double y2 = n2->y;
            double x3 = n3->x;
            double y3 = n3->y;
            double theta1 = 1; // Temperature au noeud 1
            double theta2 = 1; // Temperature au noeud 2
            double theta3 = 1; // Temperature au noeud 3
            double thetak = 0;
            // double cons = surfaceRegion[elements[i]->region][0];// alpha

            gmm::dense_matrix<double> Ke(3,3); // Second terme avec la conductivité sous forme matricielle
            gmm::dense_matrix<double> KTe(3,1); // Le nouveau terme ne dépend pas de J, on le prend comme un vecteur qu'on calculera à part
            gmm::dense_matrix<double> alpha(2,2); // matrice alpha

            alpha(0,0) = surfaceRegion[elements[i]->region][1];
            alpha(0,1) = surfaceRegion[elements[i]->region][2];
            alpha(1,0) = surfaceRegion[elements[i]->region][3];
            alpha(1,1) = surfaceRegion[elements[i]->region][4];

            gmm::dense_matrix<double> beta(2,2); // matrice beta pour la conductivité

            beta(0,0) = surfaceRegion[elements[i]->region][5];
            beta(0,1) = surfaceRegion[elements[i]->region][6];
            beta(1,0) = surfaceRegion[elements[i]->region][7];
            beta(1,1) = surfaceRegion[elements[i]->region][8];

            double J = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);

            KTe(0,0) = ((theta2 - theta1)*(alpha(0,0)*(y3-y2)*(y3-y1) + alpha(1,1)*(x2-x3)*(x1-x3) + alpha(0,1)*((x2-x3)*(y3-y1) + (y3-y2)*(x1-x3))) + (theta3-theta1)*(alpha(0,0)*(y3-y2)*(y1-y2) + alpha(1,1)*(x2-x3)*(x2-x1) + alpha(1,0)*((x2-x3)*(y1-y2) + (x2-x1)*(y3-y2))))/6;
            KTe(1,0) = (- (theta2 - theta1)*(alpha(0,0)*(y3-y1)*(y3-y1) + alpha(1,1)*(x1-x3)*(x1-x3) + 2*alpha(0,1)*(x1-x3)*(y3-y1)) - (theta3-theta1)*(alpha(0,0)*(y3-y1)*(y1-y2) + alpha(1,1)*(x1-x3)*(x2-x1) + alpha(1,0)*((x1-x3)*(y1-y2) + (x2-x1)*(y3-y1))))/6;
            KTe(2,0) = (- (theta2 - theta1)*(alpha(0,0)*(y3-y1)*(y1-y2) + alpha(1,1)*(x2-x1)*(x1-x3) + alpha(0,1)*((x1-x3)*(y1-y2) + (y3-y1)*(x2-x1))) - (theta3-theta1)*(alpha(0,0)*(y1-y2)*(y1-y2) + alpha(1,1)*(x2-x1)*(x2-x1) + 2*alpha(1,0)*(x2-x1)*(y1-y2)))/6;


            Ke(0,0) = ((alpha(0,0)*thetak + beta(0,0))*(y3-y2)*(y3-y2) + (alpha(1,1)*thetak + beta(1,1))*(x2-x3)*(x2-x3) + 2*(alpha(1,0)*thetak + beta(1,0))*(y3-y2)*(x2-x3))/2;
            Ke(1,0) = ((alpha(0,0)*thetak + beta(0,0))*(y3-y2)*(y3-y1) + (alpha(1,1)*thetak + beta(1,1))*(x2-x3)*(x1-x3) + (alpha(1,0)*thetak + beta(1,0))*((y3-y2)*(x1-x3) + (x2-x3)*(y3-y1)))/2;
            Ke(2,0) = ((alpha(0,0)*thetak + beta(0,0))*(y3-y2)*(y1-y2) + (alpha(1,1)*thetak + beta(1,1))*(x2-x3)*(x2-x1) + (alpha(1,0)*thetak + beta(1,0))*((y3-y2)*(x2-x1) + (x2-x3)*(y1-y2)))/2;
            Ke(0,1) = Ke(1,0);
            Ke(1,1) = ((alpha(0,0)*thetak + beta(0,0))*(y3-y1)*(y3-y1) + (alpha(1,1)*thetak + beta(1,1))*(x1-x3)*(x1-x3) + 2*(alpha(1,0)*thetak + beta(1,0))*(y3-y1)*(x1-x3))/2;
            Ke(2,1) = ((alpha(0,0)*thetak + beta(0,0))*(y3-y1)*(y1-y2) + (alpha(1,1)*thetak + beta(1,1))*(x1-x3)*(x2-x1) + (alpha(1,0)*thetak + beta(1,0))*((y3-y1)*(x2-x1) + (x1-x3)*(y1-y2)))/2;
            Ke(0,2) = Ke(2,0);
            Ke(1,2) = Ke(2,1);
            Ke(2,2) = ((alpha(0,0)*thetak + beta(0,0))*(y1-y2)*(y1-y2) + (alpha(1,1)*thetak + beta(1,1))*(x2-x1)*(x2-x1) + 2*(alpha(1,0)*thetak + beta(1,0))*(y1-y2)*(x2-x1))/2;


            /*Utilisation du mapping NodesCorresp.
            Si une valeur doit être ajoutée à la ligne d'un noeud de droite, celle-ci est directement ajoutée à la ligne correspondant au noeud de gauche en vis-a-vis*/
            int num1 = NodesCorresp[n1]->num-1;
            int num2 = NodesCorresp[n2]->num-1;
            int num3 = NodesCorresp[n3]->num-1;

            Tmp(num1, n1->num-1) += (Ke(0,0) + KTe(0,0))/J;
            Tmp(num1, n2->num-1) += (Ke(0,1) + KTe(0,0))/J;
            Tmp(num1, n3->num-1) += (Ke(0,2) + KTe(0,0))/J;
            Tmp(num2, n1->num-1) += (Ke(1,0) + KTe(1,0))/J;
            Tmp(num2, n2->num-1) += (Ke(1,1) + KTe(1,0))/J;
            Tmp(num2, n3->num-1) += (Ke(1,2) + KTe(1,0))/J;
            Tmp(num3, n1->num-1) += (Ke(2,0) + KTe(2,0))/J;
            Tmp(num3, n2->num-1) += (Ke(2,1) + KTe(2,0))/J;
            Tmp(num3, n3->num-1) += (Ke(2,2) + KTe(2,0))/J;
        }
    }
    gmm::copy(Tmp,KT);
}


bool End_Criterion(std::vector<double> &theta_k, std::vector<double> &delta_theta_k)
{
    double eps = 1e-6;
    for (i=0; i<theta_k.size();i++)
    {
        if(abs(delta_theta_k(i)/theta_k(i)) > eps)
            return false;
    }

    return true;
}
