# [DU-Bii](https://du-bii.github.io/accueil) - module 5 - Méthodes et outils bioinformatiques pour l'analyse des données à haut débit

## Responsables du Module
* Olivier Kirsh - Université Paris-Diderot - olivier.kirsh@univ-paris-diderot.fr
* Olivier Rué - INRA - olivier.rue@inrae.fr

## Intervenants
* Claude Thermes - Institut de Biologie Intégrative de la Cellule - Claude.THERMES@i2bc.paris-saclay.fr
* Olivier Rué - INRAE - olivier.rue@inra.fr
* Valentin Loux - INRAE - valentin.loux@inrae.fr
* Marc Deloger - Institut Gustave Roussy - marc.deloger@gustaveroussy.fr
* Nicolas Servant - Institut Curie - nicolas.servant@curie.fr
* Mélanie Petera - INRAE - melanie.petera@inrae.fr
* Binta Dieme - Université Clermont d’Auvergne (UCA) - binta.dieme@uca.fr
* Thibaut Léger - Institut Jacques Monod - thibaut.leger@ijm.fr
* Matthias Zytnicki - INRAE - matthias.zytnicki@inrae.fr

## Calendrier 2020

Calendrier du DU-Bii : <https://du-bii.github.io/accueil/img/planning_DUBii_2020.png>

## Synopsis

### [Séance 1](https://github.com/DU-Bii/module-5-Methodes-Outils/tree/master/seance1)

- Date : 9 mars 2020
- Horaires : 10h00 - 17h30
- Intervenants : Olivier Rué et Valentin Loux
- Titre : First steps with NGS data


| Supports | Formats |
|--------------------------------------------------|--------|
| INTRODUCTION AU SÉQUENÇAGE À HAUT DÉBIT POUR LA GÉNOMIQUE 1 | [[pdf](seance1/20200309_THERMES_1.pdf)] [[pptx](seance1/20200309_THERMES_1.pptx)]  |
| INTRODUCTION AU SÉQUENÇAGE À HAUT DÉBIT POUR LA GÉNOMIQUE 2 | [[pdf](seance1/20200309_THERMES_2.pdf)] [[pptx](seance1/20200309_THERMES_2.pptx)] |
| INTRODUCTION AU SÉQUENÇAGE À HAUT DÉBIT POUR LA GÉNOMIQUE 3 | [[pdf](seance1/20200309_THERMES_3.pdf)] [[pptx](seance1/20200309_THERMES_3.pptx)] |
| First steps with NGS data | [[pdf](seance1/seance1.pdf)] [[html](seance1/slides.html)] |
| Correction TP | [[html](seance1//document.html)] |

### [Séance 2](https://github.com/DU-Bii/module-5-Methodes-Outils/tree/master/seance2)

- Date : 11 mars 2020
- Horaires : 9h30 - 17h00
- Intervenants : Nicolas Servant et Marc Deloger
- Titre : RNAseq data analysis


| Supports | Formats |
|--------------------------------------------------|--------|
| State of the art | [[pdf](seance2/DUBii_State_of_the_art_of_what_can_be_done_with_RNA-seq_20200311.pdf)]
| RNAseq processing | [[pdf](seance2/processing/RNAseq_processing.pdf)] [[html](seance2/processing/RNAseq_processing.html)]  |
| what2do_with_count_table_diffan | [[pdf](seance2/R/what2do_with_count_table_diffan.pdf)] [[html](seance2/R/what2do_with_count_table_diffan.html)] |

lien vers le rapport *MultiQC* sur le serveur IFB `/shared/projects/dubii2020/data/rnaseq/rawdata/multiqc_report.html`

### [Séance 3](https://github.com/DU-Bii/module-5-Methodes-Outils/tree/master/seance3)

- Date : 4 juin  2020
- Horaires : 09h30 - 12h30 - 14h00-16h30
- Intervenants : Mélanie Pétéra, Binta Diémé.
- Titre : Metabolomics


#### Pré-requis :

##### Notions générales de métabolomique :
- Avoir suivi des vidéos introductives du [MOOC FUN Métabololique](https://www.fun-mooc.fr/courses/course-v1:cnrs+136001+session01/about  
). Ce MOCC est accessible librement (après inscription) le temps de la crise sanitaire. Suivre les vidéos S1C1 (semaine 1 cours 1) , S1C3, S1C4 et S3C6. 
Si vous ne souhaitez pas vous inscire, vous pouvez télécharger directement les vidéos à [cette adresse](https://pfem.clermont.inra.fr/pydio/public/617c6e)

Temps estimé : 40 minutes.

##### Partie pratique  :

- Pour les apprenants sous Windows, avoir installé sur leur poste les logiciel [ms-convert](http://proteowizard.sourceforge.net/download.html) et [InSilicosViewer](https://ent.uca.fr/filez/tb3n2etb)  ainsi que ce [jeu de données](https://ent.uca.fr/filez/wfm5cl4tn).

Temps estimé : 10 minutes.

- L'ensemble des TPs de cours auront lieu sous l'instance Galaxy workflow4metabolomics de [usegalaxy.fr](https://workflow4metabolomics.usegalaxy.fr)

Nous vous demandons de :
- Créer un compte sur l'instance Galaxy (menu Authentification et Enregistrement). Pour les académiques, vous devriez pouvoir utiliser le bouton "Elixir login" avec vos authentifiants institutionnels. Pour les non-académiques vous pouvez utiliser le lien "Don't have an account? Register here." sur cette même page.

- Pour ceux n'ayant jamais utilisé Galaxy, suivre le tutoriel [A short introduction to Galaxy](https://galaxyproject.github.io/training-material/topics/introduction/tutorials/galaxy-intro-short/tutorial.html) 

- Faire les premieres étapes du tutoriel [Mass spectrometry: LC-MS analysis](https://galaxyproject.github.io/training-material/topics/metabolomics/tutorials/lcms/tutorial.html)  jusqu'à la partie "1.2. Data preparation for XCMS: MSnbase readMSData" comprise. C'est ce tutoriel qui sera suivi en TP. 

Temps estimé : 30mn à 1h.

A l'issu de ces activités préparatoires, vous aurez :
- des notions de bases en métabolomique
- installé un outil de conversion de données au format propriétaire en format "open-source" qui sera utilisé lors du TP
- Créé votre compte sur l'instance Galaxy usegalaxy.fr , fait vos premiers pas sous Galaxy. 

En cas de blocage ou de question, n'hésitez pas à partager vos questions sur le channel #help de Slack.

| Supports | Formats |
|--------------------------------------------------|--------|
| Cours | [[pdf](seance3/presentation-dubii_04-06-2020.pdf)]
| TP Msconvert | [[pdf](seance3/MSconvert-dubii.pdf)]  |
| Présentation Workflow4metabolomics | [[pdf](seance3/DUBii2020_W4M.pdf)] |
| Support Galaxy Training  | [[pdf](seance3/DUBii2020_GTNsuppmat.pdf)] |



### [Séance 4](https://github.com/DU-Bii/module-5-Methodes-Outils/tree/master/seance4)

- Date : 5 juin 2020
- Horaires : 9h30 - 17h00
- Intervenants : Olivier Rué et Matthias Zytnicki
- Titre : Croisement de données


| Supports | Formats |
|--------------------------------------------------|--------|
| Slides | [[html](seance4/slides.html)]
| TP | [[html](seance4/document.html)]  |

### [Séance 5](https://github.com/DU-Bii/module-5-Methodes-Outils/tree/master/seance5)

- Date : 11 juin 2020
- Horaires : 14h00 - 17h00
- Intervenants : Valentin Loux, Cédric Midoux, Olivier Rué , hélène Chiapello
- Titre : Bonnes pratiques en bioinformatique : (essayer) d'aller vers plus de reproductibilité

#### Pré-requis :
Lors du TP nous allons utiliser le logiciel de gestion de version Git et l'interface web GitHub.
Vous devez donc absolument les avoir configuré sur **votre ordinateur** et **sur votre compte sur le cluster de l'IFB**.

Voici les étapes à suivre pour configuer cela. Nous nous reposons sur l'aide de Github pour ces étapes. Elle explique les étapes à suivre selon votre système d'exploitation (Linxu, MaxOS, Windows).



##### Sur **GitHub** :
- Créer un compte sur [GitHub](https://www.github.com) avec votre adresse professionnelle.

##### Sur votre **machine locale** :

- Download and install (Git)[https://git-scm.com/downloads] (si il n’est pas déjà installé).

Configuer les informations qui seront asociés à vos "commits" sous Git : 
- Configurer son nom d'utilisateur en local : [Setting your username in Git](https://help.github.com/en/github/using-git/setting-your-username-in-git). Indiquer votre Nom, Prénom.
- Configurer son mail : [Setting your commit email address](https://help.github.com/en/github/setting-up-and-managing-your-github-user-account/setting-your-commit-email-address). Utiliser la même adresse que celle utilisée pour créer le compte github 

Les deux étapes suivantes permettent de se connecter de façon sécurisée depuis la ligne de commande à GitHub :
- Génerer une clef ssh:[Generating a new SSH key and adding it to the ssh-agent](https://help.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent) Vous n'êtes pas obligé de mettre une "passphrase". Si vous en mettez une, bien la noter, elle vous sera demandée à chaque connexion de l'outil git à github.
- Ajouter cette clef à votre compte github :[Adding a new SSH key to your GitHub account](https://help.github.com/en/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account)

- [Tester que tout est bien configuré](https://help.github.com/en/github/authenticating-to-github/testing-your-ssh-connection)

##### Sur votre **compte ifb** :

Configurer de la même façon l’utilisation de  git et GitHub :
- Se connecter par ssh sur le cluster de l'IFB
- Configurer son nom d'utilisateur : [Setting your username in Git](https://help.github.com/en/github/using-git/setting-your-username-in-git)
-  Configurer son mail : [Setting your commit email address](https://help.github.com/en/github/setting-up-and-managing-your-github-user-account/setting-your-commit-email-address). Utiliser la même adresse que celle utilisée pour créer le compte github 
Les deux étapes suivantes permettent de se connecter de façon sécurisée depuis la ligne de commande à GitHub :
- Génerer une clef ssh:[Generating a new SSH key and adding it to the ssh-agent](https://help.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent) Vous n'êtes pas obligé de mettre une "passphrase". Si vous en mettez une, bien la noter, elle vous sera demandée à chaque connexion de l'outil git à github.
- Ajouter cette clef à votre compte github :[Adding a new SSH key to your GitHub account](https://help.github.com/en/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account)


- [Tester que tout est bien configuré](https://help.github.com/en/github/authenticating-to-github/testing-your-ssh-connection)



Vous devez, pour votre ordinateur et votre compte sur l'IFB, aller jusqu'au test de connexion à Github (commande ``ssh -T git@github.com``) qui doit être concluant :

``

ssh -T git@github.com  

Hi USER! You've successfully authenticated, but GitHub does not provide shell access.qui doit être concluant   

``

En ca de soucis, n'hésitez pas à nous solliciter sur Slack.



| Supports | Formats |
|--------------------------------------------------|--------|
| Slides ||
| TP | |


